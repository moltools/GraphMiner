#!/usr/bin/env python3

#IMPORT STATEMENTS
import collections
from collections import Counter
import csv

def combine_substr(dict_smiles:dict):
    unique_substr = []
    for len_substr in dict_smiles:
        for substr in dict_smiles[len_substr]:
            if substr not in unique_substr:
                unique_substr.append(substr)
    return unique_substr

def count_freq(all_str:list):
    count = Counter(all_str)
    return count

def perc_substr(count_substr, length:int):
    all_perc = []
    for sub in count_substr:
        csv_row = (sub, count_substr[sub], count_substr[sub]/length*100)
        all_perc.append(csv_row)
    return all_perc
    
