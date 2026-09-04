#pipeline/convergence.py
import subprocess
import numpy as np
import re

def parse_sander_out(out_file: str, key: str) -> list[float]:
    values = []
    pattern = re.compile(rf"{key}\s*=\s*(-?\d+\.\d+)")
    with open(out_file) as f:
        for line in f:
            m = pattern.search(line)
            if m:
                values.append(float(m.group(1)))
    return values
#Statistical test to determine tail convergence through standard deviation of tail / mean value of tail
def tail_stable(values: list[float], tol_pct: float = 1.0, tail_frac: float = 0.2) -> bool:
    arr = np.array(values)
    tail = arr[int(len(arr) * (1 - tail_frac)):]
    pct_stdev = 100 * tail.std() / abs(tail.mean())
    return pct_stdev < tol_pct
def strided_rmsd(prmtop: str, traj: str, ref: str, mask: str = "@CA,C,N",
                  stride: int = 10, work_dir: str = ".") -> list[float]:
    cpptraj_in = f"""
parm {prmtop}
trajin {traj} 1 last {stride}
rms ToRef {ref} {mask} out rmsd_check.dat
"""
    script_path = f"{work_dir}/rmsd_check.cpptraj"
    with open(script_path, "w") as f:
        f.write(cpptraj_in)

    subprocess.run(["cpptraj", "-i", script_path], check=True, cwd=work_dir,
                    capture_output=True, text=True)

    rmsd_values = []
    with open(f"{work_dir}/rmsd_check.dat") as f:
        next(f)  # header
        for line in f:
            rmsd_values.append(float(line.split()[1]))
    return rmsd_values

def is_rmsd_plateaued(prmtop: str, traj: str, ref: str, slope_tol: float = 0.0005,
                       stride: int = 10, work_dir: str = ".") -> bool:
    rmsd = strided_rmsd(prmtop, traj, ref, stride=stride, work_dir=work_dir)
    tail = rmsd[int(len(rmsd) * 0.5):]
    slope = np.polyfit(range(len(tail)), tail, 1)[0]
    return abs(slope) < slope_tol

def check_equilibration(equil_out, prmtop, equil_traj, ref_crd, cfg) -> bool:
    density = parse_sander_out(equil_out, "Density")
    ok_density = tail_stable(density, cfg["density_tol_pct"])
    ok_rmsd = is_rmsd_plateaued(prmtop, equil_traj, ref_crd, cfg["rmsd_slope_tol"],
                                 stride=cfg.get("rmsd_stride", 10))
    if not (ok_density and ok_rmsd):
        raise RuntimeError(
            f"Equilibration did not converge (density_ok={ok_density}, rmsd_ok={ok_rmsd}). "
            "Extend equil.in nstlim and re-run before starting production."
        )
    return True