#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem
from rdkit.Chem import AllChem
import string

def select_on_size(smile_mol:str):
    ''' 
    Only return SMILES if number of heavy atoms below 20
    
    input: 
    smile_mol - SMILES format of a molecule (str)
    
    returns:
    smile_mol - SMILES format of a molecule (str)'''
    mol1 = Chem.MolFromSmiles(smile_mol)
    heavy_atoms = mol1.GetNumHeavyAtoms()
    if heavy_atoms <= 20: #Change to 40
        return smile_mol

def combine_basic_substructures(molsmiles:str):
    replacements = {}
    return molsmiles, replacements

def return_basic_substructures(molsmile:str, repl_dicts:dict, index_dicts:dict):
    return index_dicts


# def combine_basic_substructures(molsmiles:str):
#     print(molsmiles)
#     moll = Chem.MolFromSmiles(molsmiles)
#     replacements = {}
#     # list_of_atoms = []
#     # for char in molsmiles:
#     #     if char in string.ascii_uppercase:
#     #         list_of_atoms.append(char)
#     if 'C(=O)O' in molsmiles:
#         start_index = molsmiles.index('C(=O)O')
#         atom_index = start_index
#         for char in molsmiles[:start_index]:
#             if char not in string.ascii_uppercase:
#                 atom_index -= 1
#         # PROBLEM: does not work when not last atoms.
#         replacements[atom_index] = [atom_index+1, atom_index+2]
#         patt = Chem.MolFromSmiles('C(=O)O')
#         repl = Chem.MolFromSmiles('C')
#         repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
#         mts = Chem.MolToSmiles(repl_str[0])
#     else:
#         mts = molsmiles
#     return mts, replacements

# def return_basic_substructures(molsmile:str, repl_dicts:dict, index_dicts:dict):
#     mol2 = Chem.MolFromSmiles(molsmile)
#     for subgraphlength in index_dicts:
#         for subgraph in index_dicts[subgraphlength]:
#             for value in repl_dicts:
#                 if str(value) in subgraph.split('-'):
#                     index = index_dicts[subgraphlength].index(subgraph)
#                     list_val = subgraph.split('-')
#                     sorted_list = [int(index) for index in list_val]
#                     sorted_list += (repl_dicts[value])
#                     sorted_list.sort()
#                     sorted_list_str = [str(val) for val in sorted_list]
#                     new = '-'.join(sorted_list_str)
#                     index_dicts[subgraphlength][index] = new         
#     return index_dicts
