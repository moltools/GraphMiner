#!/usr/bin/env python3

#IMPORT STATEMENTS
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
    #print(inputsmiles)
    neighbours_dict = {}
    m1 = Chem.MolFromSmiles(inputsmiles)
    for atom in range(len(nodelist)):
        neighbours_dict[nodelist[atom]] = [nodelist[x.GetIdx()] for x in m1.GetAtomWithIdx(atom).GetNeighbors()]
    return neighbours_dict


def breadth_fs(nodelist:list, neigh_dict:dict):
    '''
    Find all subgraphs that are present
    
    input:
    nodelist - list containing the indeces (int) in the string containing an atom
    neigh_dicht - dictionary containing as value (int) the index of the start atom
    and as key (list of int) the indices of the neighbouring atoms

    returns:
    subgraphs - dictionary containg as value (int) the subgraph length and as key
    (list of str) with the subgraphs as strings of indeces divided by '-'
    '''
    len_subgraph = 1
    subgraphs = {}
    subgraphs[len_subgraph] = []
    #Generate the first list of keys
    for key in neigh_dict:
        subgraphs[len_subgraph].append(str(key))
    #Add in all unqiue regions for each length of subgraphs
    for index in range(len(nodelist)):
        len_subgraph += 1
        subgraphs[len_subgraph] = []
        sorted_graphs = []
        for values in subgraphs[len_subgraph-1]:
            last_added = values.split('-')[-1] 
            for nb in neigh_dict[int(last_added)]:
                if str(nb) not in values.split('-'):
                    list_val = values.split('-')
                    list_val.append(str(nb))
                    sorted_list = list_val
                    new = '-'.join(list_val)
                    sorted_list.sort()
                    if sorted_list not in sorted_graphs:
                        subgraphs[len_subgraph].append(new)
                        sorted_graphs.append(sorted_list)
        if subgraphs[len_subgraph] == []:
            #Maybe remove even last empty one?
            break
    return subgraphs


def subgraphs_smiles(sub_graphs:dict, smilesmol:str):
    '''
    Return indicies of subgraphs back to atoms
    
    input:
    subgraphs - dictionary containg as value (int) the subgraph length and as key
    (list of str) with the subgraphs as strings of indeces divided by '-'
    smilesmol - SMILES format of a molecule (str)

    return
    mol_graphs - dictionary containg as value (int) the subgraph length and as key
    (list of str) with the subgraphs as strings of atoms
    '''
    mol_graphs = {}
    for subgraph_length in sub_graphs:
        mol_graphs[subgraph_length] = []
        for subgraphlist in sub_graphs[subgraph_length]:
            mol_subgraph = ''
            for number in subgraphlist.split('-'):
                mol_subgraph += smilesmol[int(number)]
            mol_graphs[subgraph_length].append(mol_subgraph)
    return mol_graphs