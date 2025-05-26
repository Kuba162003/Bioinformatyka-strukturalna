from Bio.PDB import MMCIFParser, Superimposer, PDBParser
from Bio.PDB.Atom import Atom
import numpy as np
import matplotlib.pyplot as plt
import math

def extract_ca_atoms(structure):
    """Wyciąga wszystkie atomy CA z danej struktury."""
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if isinstance(atom, Atom) and atom.get_name() == "CA":
                        ca_atoms.append(atom)
    return ca_atoms

def calculate_gdt_ts(fixed_atoms, moving_atoms, thresholds=[1.0]):
    """Oblicza GDT_TS pomiędzy dwoma listami atomów dla różnych progów."""
    scores = []

    for threshold in thresholds:
        matching_atoms = 0
        total_atoms = len(fixed_atoms)

        for atom1, atom2 in zip(fixed_atoms, moving_atoms):
            coord1 = atom1.get_coord()
            coord2 = atom2.get_coord()
            distance = np.linalg.norm(coord1 - coord2)
            if distance <= threshold:
                matching_atoms += 1

        gdt_ts = (matching_atoms / total_atoms) * 100
        scores.append(gdt_ts)

    avg_gdt_ts = np.mean(scores)
    return avg_gdt_ts

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
thresholds=[1.0, 2.0, 4.0, 8.0]

gdt_ts_scores = []

i = 0
for file2 in files:
    # Inicjalizuj parser
    parser = MMCIFParser()

    # Parsuj struktury
    structure1 = parser.get_structure("Structure1", file1)
    structure2 = parser.get_structure("Structure2", file2)

    # Wyciągaj atomy CA
    ca_atoms1 = extract_ca_atoms(structure1)
    ca_atoms2 = extract_ca_atoms(structure2)

    # Dopasuj atomy CA
    fixed_atoms, moving_atoms = get_matching_atoms(ca_atoms1, ca_atoms2)

    # Przeprowadź superpozycję
    super_imposer = Superimposer()
    super_imposer.set_atoms(fixed_atoms, moving_atoms)
    super_imposer.apply(structure2.get_atoms())

    # Oblicz GDT_TS
    gdt_ts = calculate_gdt_ts(fixed_atoms, moving_atoms, thresholds)
    gdt_ts_scores.append(gdt_ts)
    print(f"GDT_TS {i}: {gdt_ts:.4f}")
    i += 1

# Tworzenie wykresu słupkowego
models = [f"Model {i}" for i in range(5)]
plt.bar(models, gdt_ts_scores, color='blue')
plt.xlabel("Model")
plt.ylabel("GDT_TS (%)")
plt.title("GDT_TS dla różnych modeli")
plt.ylim(90, 100)
plt.show()
