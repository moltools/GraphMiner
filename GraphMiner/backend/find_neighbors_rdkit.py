#!/usr/bin/env python3

from rdkit import Chem

def list_nodes(insmiles:str, atom_list:list):
    '''
    Create a list of all nodes based on SMILES

    input:
    insmiles - SMILES format of a molecule (str)
    atom_list - list of all atoms (and merges) to check for in molecule]

    returns:
    node_list - list containing the indeces (int) in the string containing an atom
    '''
    node_list = []
    for index in range(len(insmiles)):
        if insmiles[index:index+2] in atom_list: 
            node_list.append(insmiles[index:index+2])
        elif insmiles[index] in atom_list:
            node_list.append(insmiles[index])
    return node_list

def rdkit_parse(inputsmiles:str, nodelist:list):
    neighbours_dict = {}
    m1 = Chem.MolFromSmiles(inputsmiles)
    for atom in range(len(nodelist)):
        neighbours_dict[nodelist[atom]] = [nodelist[x.GetIdx()] for x in m1.GetAtomWithIdx(atom).GetNeighbors()]
    return neighbours_dict






