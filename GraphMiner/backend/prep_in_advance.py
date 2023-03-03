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

# def combine_basic_substructures(molsmiles:str):
#     replacements = {}
#     return molsmiles, replacements

# def return_basic_substructures(molsmile:str, repl_dicts:dict, index_dicts:dict):
#     return index_dicts


def combine_basic_substructures(molsmiles:str):
    moll = Chem.MolFromSmiles(molsmiles)
    replacements = {}
    if 'C(=O)O' in molsmiles:
        start_index = molsmiles.index('C(=O)O')
        atom_index = start_index
        for char in molsmiles[:start_index]:
            if char not in string.ascii_uppercase:
                atom_index -= 1
        # PROBLEM: does not work when not last atoms.
        replacements[atom_index] = [atom_index+1, atom_index+2]
        patt = Chem.MolFromSmiles('C(=O)O')
        repl = Chem.MolFromSmiles('C')
        repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
        mts = Chem.MolToSmiles(repl_str[0])
    else:
        mts = molsmiles
    return mts, replacements

def return_basic_substructures(repl_dicts:dict, index_dicts:dict):
    for subgraphlength in index_dicts:
        for subgraph in index_dicts[subgraphlength]:
            for value in repl_dicts:
                if str(value) in subgraph.split('-'):
                    index = index_dicts[subgraphlength].index(subgraph)
                    list_val = subgraph.split('-')
                    sorted_list = [int(indexs) for indexs in list_val]
                    length_input = len(repl_dicts[value])
                    # print('before: ' + str(sorted_list))
                    for ind in range(len(sorted_list)):
                        if sorted_list[ind] > value:
                            sorted_list[ind] = sorted_list[ind] + length_input
                    # print('after add: ' + str(sorted_list))
                    #ADD IN OTHER INDECES AND REPLACE
                    sorted_list += (repl_dicts[value])
                    sorted_list.sort()
                    # print('after including: ' + str(sorted_list))
                    # print(' ')
                    sorted_list_str = [str(val) for val in sorted_list]
                    new = '-'.join(sorted_list_str)
                    index_dicts[subgraphlength][index] = new
    return index_dicts



    for subgraphlength in index_dicts:
        for subgraph in index_dicts[subgraphlength]:
            for value in repl_dicts:
                if str(value) in subgraph.split('-'):
                    index = index_dicts[subgraphlength].index(subgraph)
                    list_val = subgraph.split('-')
                    sorted_list = [int(indexs) for indexs in list_val]
                    for number in repl_dicts[value]:
                        if number in sorted_list:
                            for ind in range(len(sorted_list)):
                                if sorted_list[ind] >= number:
                                    sorted_list[ind] += 1         
                    sorted_list += (repl_dicts[value])
                    sorted_list.sort()
                    sorted_list_str = [str(val) for val in sorted_list]
                    new = '-'.join(sorted_list_str)
                    index_dicts[subgraphlength][index] = new
                else:
                    for number in repl_dicts[value]:
                        sorted_list = [int(index) for index in subgraph.split('-')]
                        if number in sorted_list:
                            for ind in range(len(sorted_list)):
                                if sorted_list[ind] >= number:
                                    sorted_list[ind] += 1   
                            sorted_list_str = [str(val) for val in sorted_list]
                            new = '-'.join(sorted_list_str)
                            index_dicts[subgraphlength][index] = new    