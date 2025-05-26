from itertools import combinations
import math
from Bio.PDB import PDBParser, PDBList
from Bio.PDB.Atom import Atom
import numpy as np
from typing import List


def drop_backbone(atoms:List[Atom]):
    new_list=[]
    for atom in atoms:
        if atom.name not in ['N','CA','C']:
            new_list.append(atom)
    return new_list
def drop_heteroatoms(atoms:List[Atom]):
    new_list=[]
    for atom in atoms:
        if atom.parent.id[0]==' ':
            new_list.append(atom)
    return new_list
def drop_elements(atoms:List[Atom],elements:List[str]):
    new_list=[]
    for atom in atoms:
        if atom.element not in elements:
            new_list.append(atom)
    return new_list
def is_backbone(atom1:Atom,atom2:Atom):
    backbone=['N','CA','C']
    if atom1.parent == atom2.parent:
        atoms_string=str(atom1.name)+str(atom2.name)
        if atoms_string in "".join(backbone) or atoms_string in "".join(backbone[::-1]):
            return True
    if atom1.parent.id[1]==atom2.parent.id[1]-1:
        if atom1.name=='C' and atom2.name=='N':
            return True
    if atom1.parent.id[1]==atom2.parent.id[1]+1:
        if atom1.name=='N' and atom2.name=='C':
            return True
    return False
    
def extract_all_atoms(structure):
    all_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if isinstance(atom, Atom):
                        all_atoms.append(atom)
    return all_atoms

def count_clashes(atoms:List[Atom],treshold:float|int=0.4)->int:
    clash_count:int=0

    for atom1,atom2 in combinations(atoms,2):
        if is_backbone(atom1,atom2):
            continue
        distance=atom1-atom2
        if distance<2-treshold:
            clash_count+=1
    return clash_count
def clash_count_to_score(clash_count:int,atoms:List[Atom])->float:
    return float((clash_count*1000)/len(atoms))
    
pdbl = PDBList()
fetch = pdbl.retrieve_pdb_file('4ywo', file_format='pdb')

parser = PDBParser()
structure = parser.get_structure("structure", fetch)

all_atoms = extract_all_atoms(structure)
atoms=drop_heteroatoms(all_atoms)
atoms=drop_elements(atoms,['H'])
clash_count=count_clashes(atoms)
clash_score=clash_count_to_score(clash_count,atoms)
print(f"Clash Score: {clash_score:.4f}")