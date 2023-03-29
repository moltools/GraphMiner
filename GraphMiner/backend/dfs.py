#!/usr/bin/env python3

#Import statements
import rdkit
from rdkit import Chem

def depth_fs(molecule, nodedict:dict):
    substr_list = []
    for atom in molecule.GetAtoms():
        start_atom = atom.GetAtomMapNum()
        queue = nodedict[start_atom]


    return