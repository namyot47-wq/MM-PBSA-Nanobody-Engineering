#pipeline/stages.py
import subprocess
from pathlib import Path

#This mirrors the four commands used in the amber tutorial
def run_sander(stage: str, mdin: str, prmtop: str, in_crd: str,
               ref_crd: str = None, write_traj: bool = False, work_dir: Path = Path(".")):
    out = work_dir / f"{stage}.out"
    rst = work_dir / f"{stage}.rst"
    cmd = ["sander", "-O", "-i", mdin, "-o", str(out),
           "-p", prmtop, "-c", in_crd, "-r", str(rst)]
    if write_traj:
        cmd += ["-x", str(work_dir / f"{stage}.mdcrd")]
    if ref_crd:
        cmd += ["-ref", ref_crd]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    if result.returncode != 0:
        raise RuntimeError(f"[{stage}] sander failed:\n{result.stderr[-2000:]}")
    if not out.exists() or "Total wall time" not in out.read_text():
        raise RuntimeError(f"[{stage}] did not finish cleanly, check {out}")
    return rst
#pulls files from the equilibrium_templates to run equlibritaion of complex,ligand and receptor.
def run_full_equilibration(prmtop, solvated_inpcrd, work_dir: Path):
    min_rst    = run_sander("min",     "min.in",     prmtop, solvated_inpcrd,
                             ref_crd=solvated_inpcrd, work_dir=work_dir)
    heat_rst   = run_sander("heat",    "heat1.in",   prmtop, str(min_rst),
                             ref_crd=str(min_rst), write_traj=True, work_dir=work_dir)
    density_rst = run_sander("density","density.in", prmtop, str(heat_rst),
                             ref_crd=str(heat_rst), write_traj=True, work_dir=work_dir)
    equil_rst  = run_sander("equil",   "equil.in",   prmtop, str(density_rst),
                             write_traj=True, work_dir=work_dir)
    return equil_rst