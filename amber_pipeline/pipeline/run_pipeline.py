import yaml
from pathlib import Path
from pipeline import prep, render, stages, convergence

def main(config_path="config.yaml"):
    cfg = yaml.safe_load(Path(config_path).read_text())
    work_dir = Path("work") / cfg["run_id"]
    work_dir.mkdir(parents=True, exist_ok=True)

    #Tutorial step 1: Prepping of complex, receptor and ligand
    clean_pdb = work_dir / "01_prep" / "protein_complex_clean.pdb"
    prep.strip_hetero(Path(cfg["input_pdb"]), clean_pdb, cfg["keep_residues"])
    protein_mask = prep.protein_mask_from_prmtop(str(work_dir / "protein_complex_solvated.prmtop"))

    #Tutorial step 2: From TLEAP, render the pdb inputs and equlibrate
    ctx = {**cfg, "protein_mask": protein_mask}
    for tmpl, out in [("min.in.j2", "min.in"), ("heat1.in.j2", "heat1.in"),
                       ("density.in.j2", "density.in"), ("equil.in.j2", "equil.in")]:
        render.render_input(tmpl, ctx, str(work_dir / out))

    equil_rst = stages.run_full_equilibration(
        str(work_dir / "protein_complex_solvated.prmtop"),
        str(work_dir / "protein_complex_solvated.inpcrd"),
        work_dir,
    )

    convergence.check_equilibration(
        str(work_dir / "equil.out"),
        str(work_dir / "protein_complex_solvated.prmtop"),
        str(work_dir / "equil.mdcrd"),
        str(work_dir / "protein_complex_solvated.inpcrd"),
        cfg,
    )

    render.render_input("prod.in.j2", ctx, str(work_dir / "prod.in"))
    stages.run_production(str(work_dir / "protein_complex_solvated.prmtop"), equil_rst,
                           str(work_dir / "prod.in"), work_dir,
                           n_segments=cfg.get("n_prod_segments", 1))

if __name__ == "__main__":
    main()
