#!/usr/bin/env python3

#IMPORTS
from statsmodels.stats import multitest
import pandas as pd
import numpy as np
from rdkit import Chem
from scipy.stats import hypergeom

def join_df2(list_of_df):
    '''
    Combine the separate dataframes into one and replace NaN-values

    input:
    list_of_df - list containing dataframes, the dataframes contain two columns
    - substructure in SMILES (str) and OccurrenceList[groupnumber] displaying
    in str whether the substructure occurred (1) or not (0) in a molecule at
    the index of the molecule

    returns:
    joined_df - dataframe containing three columns - substructure in SMILES
    (str) and twice an OccurrenceList[groupnumber] displaying in str
    whether the substructure occurred (1) or not (0) in a molecule at the
    index of the molecule
    '''
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    colmn = list(joined_df.columns)
    joined_df[colmn[1]] = joined_df[colmn[1]].replace(np.nan, 0)
    joined_df[colmn[2]] = joined_df[colmn[2]].replace(np.nan, 0)
    return joined_df

def hypergeometric_test_pval(total_dict: dict, input_df, grouplist):
    #M = Totaal aantal moleculen
    #s = Totaal aantal waar substructuur in zit
    #N = Moleculen in deze groep
    #k = Moleculen deze groep waar substructuur in zit
    M = sum(total_dict.values())
    group = 1
    pval = {}
    for val in total_dict.keys():
        pval[val] = []
        N = total_dict[val]
        for substr in input_df.iterrows():
            for input in substr:
                if type(input) != int:
                    k = input[group]
                    s = 0
                    for index in range(len(grouplist)):
                        s += input[index + 1]
                    pval[val].append(round(1-hypergeom.cdf(k-1,M,s,N), 2))
        group += 1
    return pval

def mul_test_corr(list_of_pval, corr_meth):
    '''
    Perform Multiple Testing Correction using input method

    input:
    list_of_pval - list containing all p-values calculated for each substructure
    in order of the substructure in the dataframe

    returns:
    tests - result of multiple testing correction as array containing [1] True/
    False, [2] corrected p-values, [3] corrected alpha for Sidak method and
    [4] corrected alpha for Multiple Testing Correction method
    '''
    tests = multitest.multipletests(pvals = list_of_pval, alpha = 0.05, method = corr_meth)
    TF_list_mtc = tests[0].tolist()
    return TF_list_mtc

def extract_signif_substr(TF_list, df_substr):
    '''
    Retrieve all substructures that are significantly differently present

    input:
    TF_list - list containing a True or False for each molecule after
    multiple testing correction
    df_substr - dataframe containing the substructure, whether it is present
    in each molecule of the group (1) or not (0), and whether the multiple
    testing correction returns True or False

    returns:
    sig_dif - a list containing all significantly differently present
    substructures in descending length
    '''
    index = 0
    sig_dif_dict = {}
    sig_dif = []
    print(df_substr)
    for TF in TF_list:
        if TF == True:
            sig_dif.append(df_substr.iloc[index][0])
            value0 = df_substr.iloc[index][1]
            value1 = df_substr.iloc[index][2]
            if value0 > value1:
                sig_dif_dict[df_substr.iloc[index][0]] = 'group 0'
            elif value1 > value0:
                sig_dif_dict[df_substr.iloc[index][0]] = 'group 1'
            if value1 == value0:
                sig_dif_dict[df_substr.iloc[index][0]] = 'same'
        index += 1
    sig_dif.sort(key=len, reverse=True)
    return sig_dif, sig_dif_dict

def create_groups_substr(dif_sig):
    new_dic = {}
    for val in dif_sig:
        placed = 'no'
        for keyval in new_dic.keys():
            if Chem.MolFromSmiles(keyval).HasSubstructMatch(Chem.MolFromSmiles(val)):
                new_dic[keyval].append(val)
                placed = 'yes'
        if placed == 'no':
            new_dic[val] = []
    return new_dic
