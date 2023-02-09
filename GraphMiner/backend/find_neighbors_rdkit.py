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

def breadth_fs(node_list:list, molsmiles:str, node_dict:dict):
    '''
    Perform Breadth First Search
    
    input:
    node_list - list containing the indeces (int) in the string containing an atom
    molsmiles - SMILES format of a molecule (str)
    node_dict - dictionary containing as key (int) the index of the atom and as value
    a list of indeces (int) of the neighboring atoms

    returns:
    vis - dictionary containing as key the index (int) of the atom and as value the 
    distance (int) to the first atom
    '''
    ###WERKT NIET AANGEZIEN HET BEIDE KANTEN OP IS IN DEZE EN IN MIJN EIGEN EEN KANT OP :(
    start_node = node_list[0]
    unvis = node_list
    curr_dist = 0
    queue = [start_node]
    distances = [0]
    vis = {}
    while len(queue) != 0:
        curr_node = queue[0]
        curr_dist = distances[0]
        vis.update({(curr_node):str(curr_dist)})
        if curr_node not in node_dict.keys():
            pass
        else:
            for node in node_dict[curr_node]:
                if node in queue:
                    continue
                elif (node) not in unvis:
                    continue
                else:
                    queue.append(node)
                    distances.append(int(curr_dist)+1)
                    unvis.remove(((node)))
        queue.pop(0)
        distances.pop(0)
    node_list.pop(0)
    for nodes in unvis:
        vis.update({nodes:'-1'})
    return vis





