#!/usr/bin/env python3

#Import statements
import rdkit
from rdkit import Chem

def depth_fs(molecule, nodedict:dict, numberslist:list):
    substr_list = []
    unvis = numberslist
    vis = []
    for atom in molecule.GetAtoms():

        start_atom = atom.GetAtomMapNum()
        queue = nodedict[start_atom]
        unvis.remove(start_atom)
        vis.append(start_atom)
        # print(unvis)
        print(queue)
        for nb in queue:
            queue.pop(0)
            print(nb)



    return