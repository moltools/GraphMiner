#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem
import string


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
            node_list.append(index)
        elif insmiles[index] in atom_list:
            node_list.append(index)
    return node_list


def rdkit_parse(inputsmiles:str, nodelist:list):
    '''
    Parsing SMILES using RDKit
    
    input:
    inputsmiles - SMILES format of a molecule (str)
    nodelist - list containing the indeces (int) in the string containing an atom

    returns:
    neighbours_dict - dictionary containing as key (int) the index of the start atom
    and as value (list of int) the indices of the neighbouring atoms
    '''
    neighbours_dict = {}
    m1 = Chem.MolFromSmiles(inputsmiles)
    for atom in range(len(nodelist)):
        neighbours_dict[nodelist[atom]] = [nodelist[x.GetIdx()] for x in m1.GetAtomWithIdx(atom).GetNeighbors()]
    return neighbours_dict
