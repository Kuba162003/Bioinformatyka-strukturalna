from Bio.PDB.PDBList import PDBList
from Bio import PDB
from Bio.PDB.vectors import Vector

def is_purine(residue):
    """Check if the residue is a purine (Adenine or Guanine)."""
    return residue.get_resname() in ["A", "G"]

def is_pyrimidine(residue):
    """Check if the residue is a pyrimidine (Cytosine, Uracil, or Thymine)."""
    return residue.get_resname() in ["C", "U", "T"]

def superimpose_rna(structure):
    # Run superimposer for RNA
    list_coords = []
    side_chain_atoms = []
    
    for model in structure:
        for chain in model:
            for res in chain:
                list_atoms = []
                side_chain = []
                
                for atom in res.get_atoms():
                    # Backbone atoms common to all RNA residues
                    if atom.get_name() in ['P', "C4'"]:
                        list_atoms.append(atom)
                    
                    # Purine-specific atoms
                    elif is_purine(res) and atom.get_name() in ['N9', 'C2', 'C6']:
                        list_atoms.append(atom)
                    
                    # Pyrimidine-specific atoms
                    elif is_pyrimidine(res) and atom.get_name() in ['N1', 'C2', 'C4']:
                        list_atoms.append(atom)
                    
                    # Other atoms, saved as "side chain" atoms
                    else:
                        side_chain.append(atom)

                list_coords.append(list_atoms)
                side_chain_atoms.append(side_chain)

                if len(list_coords) == 2:  # Only take the first 2 residues for example
                    break

    # Superimpose coordinates
    sup = PDB.Superimposer()
    sup.set_atoms(fixed=list_coords[0], moving=list_coords[1])
    print("Fixed structure coordinates:")
    print([x.get_coord() for x in list_coords[0]])
    print('------------')
    print("Moving structure coordinates:")
    print([x.get_coord() for x in list_coords[1]])

try:
    # Fetch PDB structure
    structure_name = '430D'  # Example RNA PDB code
    pdbl = PDBList()
    fetch_pdb = pdbl.retrieve_pdb_file(structure_name, file_format='pdb')

    # Parse the structure
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(structure_name, fetch_pdb)

    superimpose_rna(structure)

except FileNotFoundError:
    print("Podano zla nazwe pliku, lub istnieje problem z pobraniem pliki")
except Exception as e:
    print("Wystapil inny blad:", str(e))
