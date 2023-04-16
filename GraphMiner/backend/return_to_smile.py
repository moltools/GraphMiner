#!/usr/bin/env python3

import string
from rdkit import Chem

# def subgraphs_smiles(sub_graphs:dict, smilesmol:str):
#     '''
#     Return indicies of subgraphs back to atoms
    
#     input:
#     subgraphs - dictionary containg as key (int) the subgraph length and as value
#     (list of str) with the subgraphs as strings of indeces divided by '-'
#     smilesmol - SMILES format of a molecule (str)

#     return
#     mol_graphs - dictionary containg as key (int) the subgraph length and as value
#     (list of str) with the subgraphs as strings of atoms in SMILES format
#     '''
#     mol_graphs = {}
#     alphabet = list(string.ascii_uppercase)
#     for subgraph_length in sub_graphs:
#         mol_graphs[subgraph_length] = []
#         for subgraphlist in sub_graphs[subgraph_length]:
#             mol_subgraph = ''
#             for number in subgraphlist.split('-'):
#                 index = int(number)
#                 mol_subgraph += smilesmol[index]
#                 #Don't know if this should be in?
#                 #if index == int(subgraphlist.split('-')[-1]):
#                 #    break
#                 for index in range(index, len(smilesmol)-1):
#                     index += 1
#                     if smilesmol[index] in alphabet:
#                         break
#                     else:
#                          mol_subgraph += smilesmol[index]
#             mol_graphs[subgraph_length].append(mol_subgraph)
#     return mol_graphs

def rdkit_smiles(sub_graphs:dict, smilesmol:str):
    '''
    Convert the indices of subgraphs to SMILES subgraphs

    input:
    subgraphs - dictionary containing as key (int) the subgraph length and as value
    (list of str) with the subgraphs as strings of indeces divided by '-'
    smilesmol - SMILES format of a molecule (str)

    returns:
    mol_graphs - dictionary containing as key (int) the subgraph length and as value
    (list of str) the subgraphs as strings in SMILES format
    '''
    mol = Chem.MolFromSmiles(smilesmol)
    mol_graphs = {}
    for subgraph_length in sub_graphs:
        mol_graphs[subgraph_length] = []
        for subgraphlist in sub_graphs[subgraph_length]:
            atom_list = subgraphlist.split('-')
            int_atom_list = [int(index) for index in atom_list]
            mol_graphs[subgraph_length].append(Chem.MolFragmentToSmiles(mol, int_atom_list))
    return mol_graphs

def rdkit_smiles2(sub_graphs:dict, smilesmol):
    '''
    Convert the indices of subgraphs to SMILES subgraphs

    input:
    subgraphs - dictionary containing as key (int) the subgraph length and as value
    (set of int) with the subgraphs as integers
    smilesmol - MOL format of a molecule (mol)

    returns:
    mol_graphs - dictionary containing as key (int) the subgraph length and as value
    (list of str) the subgraphs as strings in SMILES format
    '''
    mol_graphs = {}
    atommapdict = {}
    for atom in smilesmol.GetAtoms():
        if atom.GetAtomMapNum() != atom.GetIdx():
            atommapdict[atom.GetAtomMapNum()] = atom.GetIdx()
    if len(atommapdict) >0:
        for subgraph_length in sub_graphs:
            mol_graphs[subgraph_length] = []
            for subgraphset in sub_graphs[subgraph_length]:
                for atommapnum in subgraphset:
                    if atommapnum in atommapdict.keys():
                        subgraphset.remove(atommapnum)
                        subgraphset.add(atommapdict[atommapnum])
                subgraphlist = list(subgraphset)
                mol_graphs[subgraph_length].append(
                    Chem.MolFragmentToSmiles(smilesmol, subgraphlist))
    else:
        for subgraph_length in sub_graphs:
            mol_graphs[subgraph_length] = []
            for subgraphset in sub_graphs[subgraph_length]:
                subgraphlist = list(subgraphset)
                mol_graphs[subgraph_length].append(Chem.MolFragmentToSmiles(smilesmol, subgraphlist))
    return mol_graphs
