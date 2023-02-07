#!/usr/bin/env python3

from sys import argv
import pandas as pd

def load_data():
    df = pd.read_csv('test_file.csv')
    return df

def determine_groups(dataset):
    groups = dataset.iloc[:,1]
    dif_groups = []
    for group in list(groups):
        if group in dif_groups:
            continue
        else:
            dif_groups.append(group)
    return dif_groups

def create_dict(diff_groups, data_set):
    list_of_mol = {}
    for group_num in diff_groups:
        list_of_mol[group_num] = []
        for mol in range(len(data_set)):
            if data_set.iloc[mol,1] == group_num:
                list_of_mol[group_num].append(data_set.iloc[mol,0])
    return list_of_mol


