import yaml
from pathlib import Path
from pipeline import prep, render, stages, convergence


def main(config_path="config.yaml"):
    cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))

    #1 Directory Setup:
    work_dir = Path("work") / cfg["run_id"]
    prep_dir = work_dir / "01_prep"
    equil_dir = work_dir / "02_equilibration"
    prod_dir = work_dir / "03_production"

    for d in (prep_dir, equil_dir, prod_dir):
        d.mkdir(parents=True, exist_ok=True)
    
    #Tutorial step 1: Prepping of complex, receptor and ligand
    clean_pdb = prep_dir / "protein_complex_clean.pdb"
    prep.strip_hetero(Path(cfg["input_pdb"]).resolve(), clean_pdb, cfg["keep_residues"])

    # First tleap pass: Build gas-phase complex topology, report charge
    tleap_in = prep.write_tleap_script(
        clean_pdb.name, cfg["forcefield"], cfg["water_model"], cfg["box_padding_ang"], prep_dir
    )
    tleap_stdout = prep.run_tleap(tleap_in, cwd=prep_dir)
    charge = prep.parse_tleap_charge(tleap_stdout)
    print(f"[prep] tleap-reported net charge: {charge}")
    # Second tleap pass: Dynamic water leaprc loading + solvate/neutralize
    water_rc = f"leaprc.water.{cfg['water_model'].lower()}"
    solvate_script = (
        f"source {cfg['forcefield']}\n"
        f"source {water_rc}\n\n"
        f"com = loadpdb {clean_pdb.name}\n"
        + prep.build_neutralized_solvated_script(charge, cfg["water_model"], cfg["box_padding_ang"])
    )
    solvate_in = prep_dir / "tleap_solvate.in"
    solvate_in.write_text(solvate_script, encoding="utf-8")
    prep.run_tleap(solvate_in, cwd=prep_dir)

    # Split receptor/ligand gas-phase topologies
    ranges = prep.get_chain_residue_ranges(clean_pdb)
    receptor_mask = prep.chain_mask_from_ranges(ranges, cfg["receptor_chain"])
    ligand_mask = prep.chain_mask_from_ranges(ranges, cfg["ligand_chain"])
    prep.parmed_split_complex(
        prep_dir / "protein_complex_gas.prmtop", 
        prep_dir / "protein_complex_gas.inpcrd",
        receptor_mask, ligand_mask, prep_dir,
    )

    solvated_prmtop = prep_dir / "protein_complex_solvated.prmtop"
    solvated_inpcrd = prep_dir / "protein_complex_solvated.inpcrd"
    protein_mask = prep.protein_mask_from_prmtop(str(solvated_prmtop))
    print(f"[prep] done. protein_mask = {protein_mask}")

    #Tutorial step 2: From TLEAP, render the pdb inputs and equlibrate
    ctx = {**cfg, "protein_mask": protein_mask}
    eq_templates = [
        ("min.in.j2", "min.in"), 
        ("heat.in.j2", "heat.in"),
        ("density.in.j2", "density.in"), 
        ("equil.in.j2", "equil.in")
    ]
    for tmpl, out in eq_templates:
        render.render_input(tmpl, ctx, str(equil_dir / out))

    equil_rst = stages.run_full_equilibration(
        str(solvated_prmtop), str(solvated_inpcrd), equil_dir, cfg
    )

    convergence.check_equilibration(
        str(equil_dir / "equil.out"), str(solvated_prmtop),
        str(equil_dir / "equil.mdcrd"), str(solvated_inpcrd), cfg,
    )
    print("[equil] converged.")

    if step == "equil":
        print("[pipeline] Stopping after 'equil' stage as requested.")
        return

    render.render_input("prod.in.j2", ctx, str(work_dir / "prod.in"))
    stages.run_production(
        str(solvated_prmtop), str(equil_rst), str(prod_dir / "prod.in"),
        prod_dir, cfg, n_segments=cfg.get("n_prod_segments", 1),
    )
    print("[production] done.")


