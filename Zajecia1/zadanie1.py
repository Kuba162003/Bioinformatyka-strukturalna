from math import sqrt
from Bio.PDB.PDBList import PDBList
from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
import argparse

#Dane wejsciowe
structure_name = '4YWO'
# Tworzymy parser
parser = argparse.ArgumentParser()

# Dodajemy argumenty
parser.add_argument(
    "--name", type=str, help="Structure name", required=False, default="4YWO"
)

# Parsowanie argumentów
args = parser.parse_args()

structure_name = args.name

def calculate_distance(structure):
    # Lista do przechowywania współrzędnych atomów CA
    ca_atoms = []

    # Zbieranie współrzędnych atomów CA z każdego łańcucha w strukturze
    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res.get_atoms():
                    if atom.get_name() == 'CA':
                        ca_atoms.append(atom.get_coord())

    # Inicjalizacja macierzy kontaktów wypełnionej zerami
    num_atoms = len(ca_atoms)
    contact_map = np.zeros((num_atoms, num_atoms), dtype=int)

    # Obliczanie odległości między każdym atomem CA
    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            distance = sqrt(sum((ca_atoms[i] - ca_atoms[j]) ** 2)) #Obliczanie odległości euklistycznej za pomocą wzoru 
            if distance < 8.0:
                contact_map[i][j] = 1
                contact_map[j][i] = 1 

    return contact_map

def map_plot(contact_map):
    #Rysowanie macierzy kontaktow
    plt.imshow(contact_map, cmap='binary', interpolation='nearest') 
    plt.title(f'Contact Map for {structure_name}')
    plt.xlabel('Atom Index')
    plt.ylabel('Atom Index')
    plt.show()  # Wyświetlenie wykresu

try:
    #Pobieranie struktury PDB
    pdbl = PDBList()
    fetch_pdb = pdbl.retrieve_pdb_file(structure_name, file_format='pdb')

    #Parsowanie - odczytywanie danych strukturanych bialka i zapisywanie ich do zmiennej structure
    pdb_parser = PDB.PDBParser()
    structure = pdb_parser.get_structure(structure_name, fetch_pdb)
    contact_map = calculate_distance(structure)
    map_plot(contact_map)

except FileNotFoundError:
    print("Podano zla nazwe pliku, lub istnieje problem z pobraniem pliki")
except:
    print("Wystapil inny blad")