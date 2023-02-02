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
    return dif_groups
