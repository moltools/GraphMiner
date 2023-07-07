#!/usr/bin/env python3

def depth_fs(molecule, nodedict:dict):
    '''
    Perform depth first search on a molecule to obtain all substructures

    input:
    molecule - molecule in Chem.Mol format
    nodedict - dictionary containing as key(int) the atom map numbers of the
    start atom and as value (list of int) the atom map numbers of the
    neighbouring atoms

    output:
    substr_list - list containing all substructures (set of int) present in the
    molecule
    '''
    substr_list = []
    for atom in molecule.GetAtoms():
        start_atom = atom.GetAtomMapNum()
        curr = {start_atom}
        substr_list.append(curr)
        substr_list = df_algorithm(start_atom, nodedict, substr_list, curr)
    return substr_list

def df_algorithm(atommapnum, nodedict, substr_list, curr):
    '''
    Perform depth first search on a molecule to obtain all substructures

    input:
    atommapnum - atom map number (int) of current atom
    nodedict - dictionary containing as key(int) the atom map numbers of the
    start atom and as value (list of int) the atom map numbers of the
    neighbouring atoms
    substr_list - list containing all substructures (set of int) present in the
    molecule
    curr - current substructure (set of int)

    output:
    substr_list - list containing all substructures (set of int) present in the
    molecule
    '''
    for nb in nodedict[atommapnum]:
        newcurr = curr.copy()
        newcurr.add(nb)
        if newcurr in substr_list:
            continue
        substr_list.append(newcurr)
        for part in curr:
            for nbpart in nodedict[part]:
                if nbpart not in curr:
                    NWCR = curr.copy()
                    NWCR.add(nbpart)
                    if NWCR not in substr_list:
                        substr_list = df_algorithm(nbpart, nodedict, substr_list, NWCR)
        substr_list = df_algorithm(nb, nodedict, substr_list, newcurr)
    return substr_list