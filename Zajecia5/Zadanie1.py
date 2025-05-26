from Bio.PDB import MMCIFParser, Superimposer
from Bio.PDB.Atom import Atom
import numpy as np

def extract_backbone_atoms(structure):
    """Wyciąga wszystkie atomy szkieletowe (N, CA, C) z danej struktury."""
    backbone_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if isinstance(atom, Atom) and atom.get_name() in ["N", "CA", "C"]:
                        backbone_atoms.append(atom)
    return backbone_atoms

def calculate_rmsd(fixed_atoms, moving_atoms):
    """Oblicza RMSD pomiędzy dwoma listami atomów."""
    distances = []
    for atom1, atom2 in zip(fixed_atoms, moving_atoms):
        coord1 = atom1.get_coord()
        coord2 = atom2.get_coord()
        distance = np.linalg.norm(coord1 - coord2)
        distances.append(distance ** 2)
    rmsd = np.sqrt(np.sum(distances) / len(distances))
    return rmsd

def get_matching_atoms(ref_atoms, target_atoms):
    """Dopasuj atomy na podstawie ich pozycji w sekwencji lub numeracji."""
    ref_residues = {atom.get_parent().id[1]: atom for atom in ref_atoms}
    target_residues = {atom.get_parent().id[1]: atom for atom in target_atoms}
    matching_ref_atoms = []
    matching_target_atoms = []
    for res_id in ref_residues:
        if res_id in target_residues:  # Dopasuj tylko wspólne reszty
            matching_ref_atoms.append(ref_residues[res_id])
            matching_target_atoms.append(target_residues[res_id])
    return matching_ref_atoms, matching_target_atoms

# Nazwy plików
file1 = "7rn1.cif"
files = [f"fold_2025_01_02_16_31_model_{i}.cif" for i in range(5)]

i = 0
for file2 in files:
    # Inicjalizuj parser
    parser = MMCIFParser()

    # Parsuj struktury
    structure1 = parser.get_structure("Structure1", file1)
    structure2 = parser.get_structure("Structure2", file2)

    # Wyciągaj atomy szkieletowe
    backbone_atoms1 = extract_backbone_atoms(structure1)
    backbone_atoms2 = extract_backbone_atoms(structure2)

    # Dopasuj atomy szkieletowe
    fixed_atoms, moving_atoms = get_matching_atoms(backbone_atoms1, backbone_atoms2)

    # Przeprowadź superpozycję
    super_imposer = Superimposer()
    super_imposer.set_atoms(fixed_atoms, moving_atoms)
    super_imposer.apply(structure2.get_atoms())

    # Oblicz RMSD
    rmsd = calculate_rmsd(fixed_atoms, moving_atoms)
    print(f"RMSD {i}: {rmsd:.4f}")
    i += 1
