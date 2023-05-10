#!/usr/bin/env python3

#Import statements
import rdkit
from rdkit import Chem

def depth_fs(molecule, nodedict:dict):
    substr_list = []
    # print(nodedict)
    for atom in molecule.GetAtoms():
        start_atom = atom.GetAtomMapNum()
        # print('atom', start_atom)
        curr = {start_atom}
        substr_list.append(curr)
        substr_list = df_algorithm(start_atom, nodedict, substr_list, curr)
    # print(substr_list)
    return substr_list

def df_algorithm(atommapnum, nodedict, substr_list, curr):
    for nb in nodedict[atommapnum]:
        # print('nb', nb)
        # print('curr', curr)
        newcurr = curr.copy()
        newcurr.add(nb)
        # print('newcurr', newcurr)
        if newcurr in substr_list:
            continue
        # print('hi')
        substr_list.append(newcurr)
        for part in curr:
            # print('part', part)
            for nbpart in nodedict[part]:
                # print('nbpart', nbpart)
                if nbpart not in curr:
                    # print('replaced', nbpart)
                    NWCR = curr.copy()
                    NWCR.add(nbpart)
                    if NWCR not in substr_list:
                        # print('NWCR', NWCR)
                        substr_list = df_algorithm(nbpart, nodedict, substr_list, NWCR)
        substr_list = df_algorithm(nb, nodedict, substr_list, newcurr)
    return substr_list