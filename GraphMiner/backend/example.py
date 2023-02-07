#!/usr/bin/env python3

import pandas as pd

def increment(x):
    df = pd.read_csv('test_file.csv')
    groups = df.iloc[:,1]
    dif_groups = []
    for group in list(groups):
        if group in dif_groups:
            continue
        else:
            dif_groups.append(group)
    list_of_mol = {}
    for group_num in dif_groups:
        list_of_mol[group_num] = []
        for mol in range(len(df)):
            if df.iloc[mol,1] == group_num:
                list_of_mol[group_num].append(df.iloc[mol,0])
    return list_of_mol
