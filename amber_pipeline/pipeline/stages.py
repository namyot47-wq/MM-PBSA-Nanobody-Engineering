import subprocess
from pathlib import Path

#This mirrors the four commands used in the amber tutorial
def run_md_engine(stage, mdin, prmtop, in_crd, ref_crd=None, write_traj=False,
                  work_dir=Path("."), engine="pmemd.cuda"):
    out = work_dir / f"{stage}.out"
    rst = work_dir / f"{stage}.rst"
    cmd = [engine, "-O", "-i", mdin, "-o", str(out),
           "-p", prmtop, "-c", in_crd, "-r", str(rst)]
    if write_traj:
        cmd += ["-x", str(work_dir / f"{stage}.mdcrd")]
    if ref_crd:
        cmd += ["-ref", ref_crd]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    if result.returncode != 0:
        raise RuntimeError(f"[{stage}] {engine} failed:\n{result.stderr[-2000:]}")
    if not out.exists() or "Total time" not in out.read_text():
        raise RuntimeError(f"[{stage}] did not finish cleanly, check {out}")
    return rst
#pulls files from the equilibrium_templates to run equlibritaion of complex,ligand and receptor.
def run_full_equilibration(prmtop, solvated_inpcrd, work_dir: Path):
    min_rst    = run_md_engine("min",     "min.in",     prmtop, solvated_inpcrd,
                             ref_crd=solvated_inpcrd, work_dir=work_dir, engine="sander")
    heat_rst   = run_md_engine("heat",    "heat.in",    prmtop, str(min_rst),
                             ref_crd=str(min_rst), write_traj=True, work_dir=work_dir)
    density_rst = run_md_engine("density","density.in", prmtop, str(heat_rst),
                             ref_crd=str(heat_rst), write_traj=True, work_dir=work_dir)
    equil_rst  = run_md_engine("equil",   "equil.in",   prmtop, str(density_rst),
                             write_traj=True, work_dir=work_dir)
    return equil_rst
#Final stage of the pipeline where trajectories are calculated from equlibrated complexes
def run_production(prmtop, equil_rst, mdin, work_dir, n_segments=1):
    prev_rst = equil_rst
    for i in range(n_segments):
        prev_rst = run_md_engine(f"prod_{i:03d}", mdin, prmtop, str(prev_rst),
                               write_traj=True, work_dir=work_dir)
    return prev_rst
