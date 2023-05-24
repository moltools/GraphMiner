#!/usr/bin/env python3

#IMPORT STATEMENTS
from collections import Counter
from rdkit import Chem

def combine_substr(dict_smiles:dict):
    '''
    Create a list containing all unique substructures in a molecule

    input:
    dict_smiles - dictionary containing as key (int) the subgraph length and as value
    (list of str) the subgraphs as strings in SMILES format

    returns:
    unique_substr - list containing all unique substructures as string in
    SMILES format
    '''
    unique_substr = []
    unique_mol = []
    for len_substr in dict_smiles:
        for substr in dict_smiles[len_substr]:
            if Chem.MolFromSmiles(substr) not in unique_mol:
                unique_substr.append(substr)
                unique_mol.append(Chem.MolFromSmiles(substr))
    return unique_substr

def combine_substr2(dict_smiles:dict):
    '''
    Create a list containing all unique substructures in a molecule

    input:
    dict_smiles - dictionary containing as key (int) the subgraph length and as value
    (list of str) the subgraphs as strings in SMILES format

    returns:
    unique_substr - list containing all unique substructures as string in
    SMILES format
    '''
    unique_substr = []
    for substr in dict_smiles:
        if substr not in unique_substr:
            unique_substr.append(substr)
    return unique_substr

def count_freq(all_str:list):
    '''
    count the frequency each substructure occurs in the group of molecules

    input:
    all_str - list containing all unique substructures as string in
    SMILES format

    returns:
    count - Counter() containing the substructures and the number of times
    it occurs
    '''
    count = Counter(all_str)
    return count

def perc_substr(count_substr, length:int):
    '''
    Create the rows for in the csv file with substructure, frequency and percentage

    input:
    count_substr - Counter() containing the substructures and the number of times
    it occurs
    length - number of molecules in the group

    returns:
    all_perc - list containing tuples with the substructure, frequency and percentage
    '''
    all_perc = []
    for sub in count_substr:
        csv_row = (sub, count_substr[sub], count_substr[sub]/length*100)
        all_perc.append(csv_row)
    return all_perc

def list_maker2(count_sub):
    '''
    Create rows for csv file with substructure and frequency

    input:
    dict_sub - dictionary with as keys the molecule numbers and as value a list
    of all substructures present
    count_sub - Counter() containing the substructures and the number of times
    it occurs

    returns:
    all_rows - list containing tuples with substructure and frequency
    '''
    all_rows = []
    for sub in count_sub:
        csv_row = (sub, count_sub[sub])
        all_rows.append(csv_row)
    return all_rows
    
