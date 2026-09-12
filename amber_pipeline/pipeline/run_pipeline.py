import argparse
import sys
import yaml
from pathlib import Path
from pipeline import prep, render, stages, convergence


def build_dirs(cfg):
    work_dir = Path("work") / cfg["run_id"]
    prep_dir = (work_dir / "01_prep").resolve()
    equil_dir = (work_dir / "02_equilibration").resolve()
    prod_dir = (work_dir / "03_production").resolve()
    for d in (prep_dir, equil_dir, prod_dir):
        d.mkdir(parents=True, exist_ok=True)
    return work_dir.resolve(), prep_dir, equil_dir, prod_dir


def run_prep(cfg, prep_dir):
    clean_pdb = prep_dir / "protein_complex_clean.pdb"
    prep.strip_hetero(Path(cfg["input_pdb"]).resolve(), clean_pdb, cfg["keep_residues"])

    tleap_in = prep.write_tleap_script(
        clean_pdb.name, cfg["forcefield"], cfg["water_model"], cfg["box_padding_ang"], prep_dir
    )
    tleap_stdout = prep.run_tleap(tleap_in, prep_dir)
    charge = prep.parse_tleap_charge(tleap_stdout)
    print(f"[prep] tleap-reported net charge: {charge}")

    water_rc = f"leaprc.water.{cfg['water_model'].lower()}"
    solvate_script = (
        f"source {cfg['forcefield']}\n"
        f"source {water_rc}\n\n"
        f"com = loadpdb {clean_pdb.name}\n"
        + prep.build_neutralized_solvated_script(charge, cfg["water_model"], cfg["box_padding_ang"])
        + "quit\n"
    )
    solvate_in = prep_dir / "tleap_solvate.in"
    solvate_in.write_text(solvate_script, encoding="utf-8")
    prep.run_tleap(solvate_in, prep_dir)

    ranges = prep.get_chain_residue_ranges(clean_pdb)
    receptor_mask = prep.chain_mask_from_ranges(ranges, cfg["receptor_chain"])
    ligand_mask = prep.chain_mask_from_ranges(ranges, cfg["ligand_chain"])
    prep.parmed_split_complex(
        prep_dir / "protein_complex_gas.prmtop",
        prep_dir / "protein_complex_gas.inpcrd",
        receptor_mask, ligand_mask, prep_dir,
    )

    protein_mask = prep.protein_mask_from_prmtop(str(prep_dir / "protein_complex_solvated.prmtop"))
    print(f"[prep] done. protein_mask = {protein_mask}")


def run_equil(cfg, prep_dir, equil_dir):
    solvated_prmtop = prep_dir / "protein_complex_solvated.prmtop"
    solvated_inpcrd = prep_dir / "protein_complex_solvated.inpcrd"
    if not solvated_prmtop.exists():
        sys.exit(f"[equil] missing {solvated_prmtop} — run --step prep first")

    protein_mask = prep.protein_mask_from_prmtop(str(prep_dir / "protein_complex_solvated.prmtop"))
    ctx = {**cfg, "protein_mask": protein_mask}
    eq_templates = [
        ("min.in.j2", "min.in"),
        ("heat.in.j2", "heat.in"),
        ("density.in.j2", "density.in"),
        ("equil.in.j2", "equil.in"),
    ]
    for tmpl, out in eq_templates:
        render.render_input(tmpl, ctx, str(equil_dir / out))

    stages.run_full_equilibration(str(solvated_prmtop), str(solvated_inpcrd), equil_dir)

    convergence.check_equilibration(
        str(equil_dir / "equil.out"), str(solvated_prmtop),
        str(equil_dir / "equil.mdcrd"), str(solvated_inpcrd), cfg,
    )
    print("[equil] converged.")


def run_production(cfg, prep_dir, equil_dir, prod_dir):
    solvated_prmtop = prep_dir / "protein_complex_solvated.prmtop"
    equil_rst = equil_dir / "equil.rst"
    if not equil_rst.exists():
        sys.exit(f"[production] missing {equil_rst} — run --step equil first")

    protein_mask = prep.protein_mask_from_prmtop(str(solvated_prmtop))
    ctx = {**cfg, "protein_mask": protein_mask}
    render.render_input("prod.in.j2", ctx, str(prod_dir / "prod.in"))
    stages.run_production(
        str(solvated_prmtop), str(equil_rst), str(prod_dir / "prod.in"),
        prod_dir, n_segments=cfg.get("n_prod_segments", 1),
    )
    print("[production] done.")


def main():
    parser = argparse.ArgumentParser(description="AMBER MM-PBSA pipeline")
    parser.add_argument("config_path", nargs="?", default="config.yaml")
    parser.add_argument("--step", choices=["prep", "equil", "production", "all"], default="all")
    args = parser.parse_args()

    cfg = yaml.safe_load(Path(args.config_path).read_text(encoding="utf-8"))
    _work_dir, prep_dir, equil_dir, prod_dir = build_dirs(cfg)

    if args.step in ("prep", "all"):
        run_prep(cfg, prep_dir)
    if args.step in ("equil", "all"):
        run_equil(cfg, prep_dir, equil_dir)
    if args.step in ("production", "all"):
        run_production(cfg, prep_dir, equil_dir, prod_dir)


if __name__ == "__main__":
    main()

