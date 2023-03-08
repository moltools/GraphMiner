#!/usr/bin/env python3

from sys import argv
import pandas as pd

def load_data(input_num:int):
    '''
    Load in the csv file

    returns: 
    df - pandas dataframe of the csv file
    '''
    df = pd.read_csv(argv[input_num], sep = ';')
    return df

def determine_groups(dataset):
    '''
    Determine which groups are present in the dataset
    
    input:
    dataset - pandas dataframe containing SMILES in first column and group number in second column

    returns:
    dif_groups - list of the group numbers that are present in the dataset
    '''
    groups = dataset.iloc[:,1]
    dif_groups = []
    for group in list(groups):
        if group in dif_groups:
            continue
        else:
            dif_groups.append(group)
    return dif_groups

def create_dict(diff_groups:list, data_set):
    '''
    Create a dictionairy of the input data, to separate SMILES on group
    
    input:
    diff_groups - list of the group numbers that are present in the data
    data_set - pandas dataframe containing SMILES in first column and group number in second column
    
    returns:
    list_of_mol - dictionary with as key the group number and as value a list of all SMILES in that group'''
    list_of_mol = {}
    for group_num in diff_groups:
        list_of_mol[group_num] = []
        for mol in range(len(data_set)):
            if data_set.iloc[mol,1] == group_num:
                list_of_mol[group_num].append(data_set.iloc[mol,0])
    return list_of_mol


