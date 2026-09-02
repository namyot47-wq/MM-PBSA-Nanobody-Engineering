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