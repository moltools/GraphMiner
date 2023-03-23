#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem


# def list_nodes(insmiles:str, atom_list:list):
#     '''
#     Create a list of all nodes based on SMILES

#     input:
#     insmiles - SMILES format of a molecule (str)
#     atom_list - list of all atoms (and merges) to check for in molecule]

#     returns:
#     node_list - list containing the indeces (int) in the string containing an atom
#     '''
#     node_list = []
#     for index in range(len(insmiles)):
#         if insmiles[index:index+2] in atom_list: 
#             node_list.append(index)
#         elif insmiles[index] in atom_list:
#             node_list.append(index)
#     return node_list


def rdkit_parse(inputsmiles:str):
    '''
    Parsing SMILES using RDKit
    
    input:
    inputsmiles - SMILES format of a molecule (str)
    nodelist - list containing the indeces (int) in the string containing an atom

    returns:
    neighbours_dict - dictionary containing as key (int) the index of the start atom
    and as value (list of int) the indices of the neighbouring atoms
    '''
    print(inputsmiles)
    neighbours_dict = {}
    m1 = Chem.MolFromSmiles(inputsmiles)
    num_heavy_atoms = m1.GetNumHeavyAtoms()
    heavy_atoms_list = list(range(0,num_heavy_atoms))
    for atom in range(num_heavy_atoms):
        neighbours_dict[heavy_atoms_list[atom]] = [heavy_atoms_list[x.GetIdx()] for x in m1.GetAtomWithIdx(atom).GetNeighbors()]
    print(neighbours_dict)
    return neighbours_dict, heavy_atoms_list


def rdkit_parse_atommap(inputmol: str):
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
    num_heavy_atoms = inputmol.GetNumHeavyAtoms()
    heavy_atoms_list = list(range(0, num_heavy_atoms))
    for atom in range(num_heavy_atoms):
        print('ATOM')
        print(atom.GetAtomMapNum())
        neighbours_dict[heavy_atoms_list[atom]] = [heavy_atoms_list[x.GetIdx()]
                                                   for x in inputmol.GetAtomWithIdx(
                atom).GetNeighbors()]
    print(neighbours_dict)
    return neighbours_dict, heavy_atoms_list