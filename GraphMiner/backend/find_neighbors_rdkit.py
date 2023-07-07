#!/usr/bin/env python3

def rdkit_parse_atommap(inputmol: str):
    '''
    Obtaining neighbors of atoms using RDKit and AtomMapNumbers

    input:
    inputmol - SMILES format of a molecule (str)

    returns:
    neighbours_dict - dictionary containing as key (int) the AtomMapNumber of
    the start atom and as value (list of int) the AtomMapNumbers of the
    neighbouring atoms
    atommapnumlist - list containing all atommapnumbers of the atoms in the
    molecule
    '''
    neighbours_dict = {}
    atommapnumlist = []
    for atom in inputmol.GetAtoms():
        atommapnumlist.append(atom.GetAtomMapNum())
        neighbours_dict[atom.GetAtomMapNum()] = \
            [x.GetAtomMapNum() for x in atom.GetNeighbors()]
    return neighbours_dict, atommapnumlist
