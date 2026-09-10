# pipeline/prep.py
import re
import subprocess
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select

def strip_hetero(input_pdb: Path, output_pdb: Path, keep: list[str] = None):
    keep = keep or []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(input_pdb))

    class ProteinOnly(Select):
        def accept_residue(self, residue):
            if residue.id[0] == " ":
                return True
            return residue.resname in keep
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_pdb), ProteinOnly())     

def write_tleap_script(complex_pdb, forcefield, water_model, box_padding, out_dir: Path):
    script = f"""
source {forcefield}
source leaprc.water.tip3p
 
com = loadpdb {complex_pdb}
 
set default PBRadii mbondi2
 
saveamberparm com protein_complex_gas.prmtop protein_complex_gas.inpcrd
 
charge com
quit
"""
    script_path = out_dir / "tleap.in"
    script_path.write_text(script)
    return script_path

def run_tleap(script_path: Path, work_dir: Path) -> str:
    result = subprocess.run(
        ["tleap", "-f", str(script_path.name)],
        cwd=str(work_dir),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"tleap failed:\n{result.stdout}\n{result.stderr}")
    return result.stdout

def parse_tleap_charge(tleap_stdout: str) -> float:
    match = re.search(r"Total (?:unperturbed )?charge:\s*(-?\d+\.\d+)", tleap_stdout)
    if not match:
        raise ValueError("Could not find total charge in tleap output")
    return float(match.group(1))
 
 
def build_neutralized_solvated_script(complex_prmtop_charge: float, water_model, padding):
    ion_lines = ""
    if abs(complex_prmtop_charge) > 1e-3:
        ion = "Cl-" if complex_prmtop_charge > 0 else "Na+"
        n_ions = round(abs(complex_prmtop_charge))
        ion_lines = f"addIons2 com {ion} {n_ions}\n"
    return f"""
{ion_lines}
solvatebox com {water_model} {padding}
saveamberparm com protein_complex_solvated.prmtop protein_complex_solvated.inpcrd
"""
 
 
def build_solvated_system(clean_pdb: Path, cfg: dict, work_dir: Path) -> Path:
    clean_pdb = clean_pdb.resolve()
    # Step 1: build gas-phase system to determine net charge
    gas_script = write_tleap_script(
        complex_pdb=str(clean_pdb),
        forcefield=cfg["forcefield"],
        water_model=cfg["water_model"],
        box_padding=cfg["box_padding_ang"],
        out_dir=work_dir,
    )
    gas_stdout = run_tleap(gas_script, work_dir)
    charge = parse_tleap_charge(gas_stdout)
 
    # Step 2: build solvated + neutralized system using that charge
    solvate_block = build_neutralized_solvated_script(
        complex_prmtop_charge=charge,
        water_model=cfg["water_model"],
        padding=cfg["box_padding_ang"],
    )
    solvate_script_path = work_dir / "tleap_solvate.in"
    solvate_script_path.write_text(f"""
source {cfg["forcefield"]}
source leaprc.water.tip3p
 
com = loadpdb {clean_pdb}
set default PBRadii mbondi2
{solvate_block}
quit
""")
    run_tleap(solvate_script_path, work_dir)
 
    return work_dir / "protein_complex_solvated.prmtop"
 

def parmed_split_complex(complex_prmtop, complex_inpcrd, receptor_mask, ligand_mask, out_dir):
    import parmed as pmd

    complex_parm = pmd.load_file(str(complex_prmtop), str(complex_inpcrd))

    receptor = complex_parm[receptor_mask]   
    ligand   = complex_parm[ligand_mask]      

    receptor.save(str(out_dir / "receptor_gas.prmtop"), overwrite=True)
    receptor.save(str(out_dir / "receptor_gas.inpcrd"), overwrite=True)
    ligand.save(str(out_dir / "ligand_gas.prmtop"), overwrite=True)
    ligand.save(str(out_dir / "ligand_gas.inpcrd"), overwrite=True)

def protein_mask_from_prmtop(prmtop_path: str) -> str:
    import parmed as pmd
    parm = pmd.load_file(prmtop_path)
    protein_residues = [r.idx + 1 for r in parm.residues if r.name not in ("WAT", "Na+", "Cl-")]
    return f":1-{max(protein_residues)}"
