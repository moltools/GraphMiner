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
    '''
    Replace acid group by C and store replacements

    input:
    molsmiles - SMILES format of a molecule (str)

    returns:
    mts - SMILES format of a molecule, with combined substructures (str)
    replacements - dictionary containing as key (int) the index of the remaining
    atom of the substructure and as value (list of int) the indeces of the replaced atoms
    '''
    moll = Chem.MolFromSmiles(molsmiles)
    mts = molsmiles
    replacements = {}
    list_of_C = [Cidx[0] for Cidx in moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
    if moll.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
        idx_COO = moll.GetSubstructMatches(Chem.MolFromSmiles('C(=O)O'))
        for idx__COO in idx_COO:
            new_list = []
            for indiv_idx in idx__COO:
                if indiv_idx in list_of_C:
                    start_idx = indiv_idx
                elif indiv_idx not in list_of_C:
                    new_list.append(indiv_idx)
        replacements[start_idx] = new_list
        #Actual Replacement
        patt = Chem.MolFromSmiles('C(=O)O')
        repl = Chem.MolFromSmiles('C')
        repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
        #ONLY WORKS FOR ONE CHAR PER MOLECULE
        mts = Chem.MolToSmiles(repl_str[0])
    ###UGH WERKT NOG NIETT
    if moll.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
        idx_CO = moll.GetSubstructMatches(Chem.MolFromSmiles('CO'))
        for idx__CO in idx_CO:
            new_list = []
            for indiv_idx in idx__CO:
                if indiv_idx in replacements.keys():
                    continue
                if indiv_idx in list_of_C:
                    start_idx = indiv_idx
                elif indiv_idx not in list_of_C:
                    new_list.append(indiv_idx)
        replacements[start_idx] = new_list
        mts = molsmiles
    print(replacements)
    return mts, replacements


def return_basic_substructures(repl_dicts:dict, index_dicts:dict):
    '''
    Add in substrucures again and

    input:
    repl_dicts - dictionary containing as key (int) the index of the remaining
    atom of the substructure and as value (list of int) the indeces of the replaced atoms
    index_dicts - dictionary containg as key (int) the subgraph length and as value
    (list of str) with the subgraphs as strings of indeces divided by '-'
    
    returns:
    index_dicts - dictionary containg as key (int) the subgraph length and as value
    (list of str) with the subgraphs as strings of indeces divided by '-', but with 
    substructures indeces placed back in again
    '''
    for subgraphlength in index_dicts:
        for subgraph in index_dicts[subgraphlength]:
            for value in repl_dicts:
                length_input = len(repl_dicts[value])
                index = index_dicts[subgraphlength].index(subgraph)
                sorted_list = [int(idx) for idx in subgraph.split('-')]
                #Replace all values that are higher than the one for which the 
                #subgraph will be added in
                for ind in range(len(sorted_list)):
                    if sorted_list[ind] > value:
                        sorted_list[ind] = sorted_list[ind] + length_input
                #Add in the earlier combined subgraphs
                if str(value) in subgraph.split('-'):
                    sorted_list += (repl_dicts[value])
                    sorted_list.sort()
                sorted_list_str = [str(val) for val in sorted_list]
                new = '-'.join(sorted_list_str)
                index_dicts[subgraphlength][index] = new
    return index_dicts