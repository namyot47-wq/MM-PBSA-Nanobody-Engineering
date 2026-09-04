# pipeline/prep.py
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

receptor = com 
ligand = com

set default PBRadii mbondi2

saveamberparm com {out_dir}/protein_complex_gas.prmtop {out_dir}/protein_complex_gas.inpcrd

charge com
"""
    (out_dir / "tleap.in").write_text(script)
    return out_dir / "tleap.in"

def parmed_split_complex(complex_prmtop, complex_inpcrd, receptor_mask, ligand_mask, out_dir):
    import parmed as pmd

    complex_parm = pmd.load_file(str(complex_prmtop), str(complex_inpcrd))

    receptor = complex_parm[receptor_mask]   
    ligand   = complex_parm[ligand_mask]      

    receptor.save(str(out_dir / "receptor_gas.prmtop"), overwrite=True)
    receptor.save(str(out_dir / "receptor_gas.inpcrd"), overwrite=True)
    ligand.save(str(out_dir / "ligand_gas.prmtop"), overwrite=True)
    ligand.save(str(out_dir / "ligand_gas.inpcrd"), overwrite=True)


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
def protein_mask_from_prmtop(prmtop_path: str) -> str:
    import parmed as pmd
    parm = pmd.load_file(prmtop_path)
    protein_residues = [r.idx + 1 for r in parm.residues if r.name not in ("WAT", "Na+", "Cl-")]
    return f":1-{max(protein_residues)}"
