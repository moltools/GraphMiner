#!/usr/bin/env python3

#Import statements
import rdkit
from rdkit import Chem

def depth_fs(molecule, nodedict:dict, numberslist:list):
    ###WERKT NOG NIET VOOR CYCLO's
    substr_list = []
    nodedict2 = {}
    for atom in molecule.GetAtoms():
        nodedict2[atom.GetAtomMapNum()] = []
        for nb in nodedict[atom.GetAtomMapNum()]:
            if nb > atom.GetAtomMapNum():
                nodedict2[atom.GetAtomMapNum()] += [nb]
    for atom in molecule.GetAtoms():
        start_atom = atom.GetAtomMapNum()
        curr = {start_atom}
        substr_list.append(curr)
        substr_list = df_algorithm(start_atom, nodedict2, substr_list, curr)
    return substr_list

def df_algorithm(atommapnum, nodedict, substr_list, curr):
    for nb in nodedict[atommapnum]:
        newcurr = curr.copy()
        newcurr.add(nb)
        substr_list.append(newcurr)
        substr_list = df_algorithm(nb, nodedict, substr_list, newcurr)
    return substr_list