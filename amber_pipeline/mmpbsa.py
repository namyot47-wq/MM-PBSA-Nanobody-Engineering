#!/usr/bin/env python3
# dG_bind = <G_complex> - <G_receptor> - <G_ligand>, single-trajectory MM-GBSA.
# Requires AmberTools (cpptraj, sander) on PATH.

import argparse
import configparser
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--complex-prmtop", required=True, type=Path)
    p.add_argument("--receptor-prmtop", required=True, type=Path)
    p.add_argument("--ligand-prmtop", required=True, type=Path)
    p.add_argument("--trajectory", required=True, type=Path)
    p.add_argument("--receptor-mask", required=True)
    p.add_argument("--ligand-mask", required=True)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--outdir", default=Path("mmpbsa_run"), type=Path)
    p.add_argument("--cpptraj", default=shutil.which("cpptraj") or "cpptraj")
    p.add_argument("--sander", default=shutil.which("sander") or "sander")
    return p.parse_args()


def load_config(path: Path) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    cfg.read(path)
    required = {
        "general": ("start_frame", "end_frame", "interval"),
        "gb": ("igb", "saltcon"),
    }
    for section, keys in required.items():
        if section not in cfg:
            raise ValueError(f"Config missing [{section}] section")
        for k in keys:
            if k not in cfg[section]:
                raise ValueError(f"Config [{section}] missing key '{k}'")
    return cfg


def run(cmd, cwd=None):
    print(f"  $ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise RuntimeError(f"Backend command failed: {cmd[0]}")
    return result


# --- backend: cpptraj -------------------------------------------------

def strip_trajectory(cpptraj_bin, complex_prmtop, trajectory,
                      mask, start, stop, interval, out_traj, workdir):
    script = workdir / f"cpptraj_{out_traj.stem}.in"
    script.write_text(
        f"parm {complex_prmtop}\n"
        f"trajin {trajectory} {start} {stop} {interval}\n"
        f"strip !({mask})\n"
        f"trajout {out_traj}\n"
        f"go\n"
    )
    run([cpptraj_bin, "-i", str(script)], cwd=workdir)
    return out_traj


# --- backend: sander ----------------------------------------------------

def write_sander_gb_input(igb, saltcon, path):
    path.write_text(
        "Single-point GB energy for each frame of a trajectory\n"
        "&cntrl\n"
        " imin=5, ntx=1, irest=0,\n"
        f" igb={igb}, saltcon={saltcon},\n"
        " ntb=0, cut=999.0,\n"
        "/\n"
    )
    return path


def run_sander_energy(sander_bin, prmtop, stripped_traj, mdin, mdout, workdir):
    run([
        sander_bin, "-O",
        "-i", str(mdin),
        "-p", str(prmtop),
        "-c", str(stripped_traj),
        "-y", str(stripped_traj),
        "-o", str(mdout),
    ], cwd=workdir)
    return mdout


# --- frontend: parse + aggregate ----------------------------------------

ENERGY_TERMS = ("VDWAALS", "EEL", "EGB", "ESURF")
_TERM_RE = {t: re.compile(rf"\b{t}\s*=\s*(-?\d+\.\d+)") for t in ENERGY_TERMS}


def parse_mdout(mdout_path: Path):
    text = mdout_path.read_text()
    values = {t: [float(m) for m in p.findall(text)] for t, p in _TERM_RE.items()}
    n_frames = len(values["EEL"])
    if n_frames == 0:
        raise RuntimeError(f"No energy frames parsed from {mdout_path}")
    return values, n_frames


def summarize(values):
    summary = {}
    for term, series in values.items():
        mean = statistics.mean(series)
        sd = statistics.stdev(series) if len(series) > 1 else 0.0
        sem = sd / (len(series) ** 0.5) if series else 0.0
        summary[term] = {"mean": mean, "stdev": sd, "sem": sem, "series": series}
    total_series = [sum(vals) for vals in zip(*values.values())]
    summary["TOTAL"] = {
        "mean": statistics.mean(total_series),
        "stdev": statistics.stdev(total_series) if len(total_series) > 1 else 0.0,
        "sem": (statistics.stdev(total_series) / (len(total_series) ** 0.5)
                if len(total_series) > 1 else 0.0),
        "series": total_series,
    }
    return summary


def main():
    args = parse_args()
    cfg = load_config(args.config)
    args.outdir.mkdir(parents=True, exist_ok=True)

    start = cfg.getint("general", "start_frame")
    stop = cfg.getint("general", "end_frame")
    interval = cfg.getint("general", "interval")
    igb = cfg.get("gb", "igb")
    saltcon = cfg.get("gb", "saltcon")
    complex_mask = cfg.get("general", "complex_mask", fallback="!(:WAT,Na+,Cl-,K+)")

    systems = {
        "complex": (args.complex_prmtop, complex_mask),
        "receptor": (args.receptor_prmtop, args.receptor_mask),
        "ligand": (args.ligand_prmtop, args.ligand_mask),
    }

    results = {}
    for name, (prmtop, mask) in systems.items():
        stripped_traj = args.outdir / f"{name}_stripped.nc"
        strip_trajectory(args.cpptraj, args.complex_prmtop, args.trajectory,
                          mask, start, stop, interval, stripped_traj, args.outdir)

        mdin = write_sander_gb_input(igb, saltcon, args.outdir / f"{name}_gb.in")
        mdout = run_sander_energy(args.sander, prmtop, stripped_traj, mdin,
                                   args.outdir / f"{name}.mdout", args.outdir)

        values, n_frames = parse_mdout(mdout)
        results[name] = summarize(values)
        print(f"[{name}] parsed {n_frames} frames")

    dG_series = [
        c - r - l for c, r, l in zip(
            results["complex"]["TOTAL"]["series"],
            results["receptor"]["TOTAL"]["series"],
            results["ligand"]["TOTAL"]["series"],
        )
    ]
    dG_mean = statistics.mean(dG_series)
    dG_sd = statistics.stdev(dG_series) if len(dG_series) > 1 else 0.0
    dG_sem = dG_sd / (len(dG_series) ** 0.5) if dG_series else 0.0

    report_path = args.outdir / "FINAL_RESULTS_simple_mmpbsa.dat"
    with open(report_path, "w") as f:
        f.write("Simple single-trajectory MM-GBSA results (no entropy term)\n")
        f.write("All units kcal/mol.\n\n")
        for name in ("complex", "receptor", "ligand"):
            f.write(f"{name.upper()}:\n")
            for term in (*ENERGY_TERMS, "TOTAL"):
                s = results[name][term]
                f.write(f"  {term:10s} {s['mean']:12.4f}  "
                        f"(stdev {s['stdev']:8.4f}, sem {s['sem']:7.4f})\n")
            f.write("\n")
        f.write(f"DELTA G binding = {dG_mean:.4f} +/- {dG_sd:.4f} "
                f"(sem {dG_sem:.4f}) kcal/mol\n")

    print(f"\nWrote {report_path}")
    print(f"DELTA G binding = {dG_mean:.4f} +/- {dG_sd:.4f} kcal/mol (sem {dG_sem:.4f})")


if __name__ == "__main__":
    main()