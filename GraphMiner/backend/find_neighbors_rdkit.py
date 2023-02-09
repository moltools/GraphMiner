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
    # for index in range(len(insmiles)):
    #     if insmiles[index:index+2] in atom_list: 
    #         node_list.append(insmiles[index:index+2])
    #     elif insmiles[index] in atom_list:
    #         node_list.append(insmiles[index])
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
    neighbours_dict - dictionary containing as value (int) the index of the start atom
    and as key (list of int) the indices of the neighbouring atoms
    '''
    print(inputsmiles)
    neighbours_dict = {}
    m1 = Chem.MolFromSmiles(inputsmiles)
    for atom in range(len(nodelist)):
        neighbours_dict[nodelist[atom]] = [nodelist[x.GetIdx()] for x in m1.GetAtomWithIdx(atom).GetNeighbors()]
    return neighbours_dict

def breadth_fs(nodelist:list, smilesmol:str, neigh_dict:dict):
    len_subgraph = 1
    subgraphs = {}
    subgraphs[len_subgraph] = []
    for key in neigh_dict:
        subgraphs[len_subgraph].append(key)
    for index in range(len(nodelist)):
        len_subgraph += 1
        subgraphs[len_subgraph] = []
        for values in subgraphs[len_subgraph-1]:
            print(str(values)[-1])
        
        
    print(subgraphs)
    return




