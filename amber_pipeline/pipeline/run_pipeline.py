import argparse
import yaml
from pathlib import Path
from pipeline import prep, render, stages, convergence


def main():
    args = parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text())
    run_id = args.run_id or cfg["run_id"]
    work_dir = Path("work") / run_id
    (work_dir / "01_prep").mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    #Tutorial step 1: Prepping of complex, receptor and ligand
    clean_pdb = work_dir / "01_prep" / "protein_complex_clean.pdb"
    clean_pdb.parent.mkdir(parents=True, exist_ok=True)
    prep.strip_hetero(Path(cfg["input_pdb"]), clean_pdb, cfg["keep_residues"])
    solvated_prmtop = prep.build_solvated_system(clean_pdb, cfg, work_dir)
    protein_mask = prep.protein_mask_from_prmtop(str(solvated_prmtop))

    if step == "prep":
        print(f"Prep complete. Solvated system written to: {solvated_prmtop}")
        return

    #Tutorial step 2: From TLEAP, render the pdb inputs and equlibrate
    ctx = {**cfg, "protein_mask": protein_mask}
    for tmpl, out in [("min.in.j2", "min.in"), ("heat.in.j2", "heat1.in"),
                       ("density.in.j2", "density.in"), ("equil.in.j2", "equil.in")]:
        render.render_input(tmpl, ctx, str(work_dir / out))

    equil_rst = stages.run_full_equilibration(
        str(work_dir / "protein_complex_solvated.prmtop"),
        str(work_dir / "protein_complex_solvated.inpcrd"),
        work_dir,
    )

    if step == "equil":
        print(f"Equilibration complete. Restart file: {equil_rst}")
        return

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
    parser = argparse.ArgumentParser(description="Run the MM-PBSA nanobody prep/equil/prod pipeline")
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument(
        "--step",
        choices=["prep", "equil", "all"],
        default="all",
        help="Which stage to run through: 'prep' (clean + solvate only), "
             "'equil' (also run equilibration), or 'all' (full pipeline through production)",
    )
    args = parser.parse_args()
    main(config_path=args.config, step=args.step)
