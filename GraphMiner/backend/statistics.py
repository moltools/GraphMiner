#!/usr/bin/env python3

#IMPORTS
from statsmodels.stats import multitest
import pandas as pd
import numpy as np
from rdkit import Chem
from scipy.stats import hypergeom

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
                print(input)
                if type(input) != int:
                    k = input[group+1]
                    s = 0
                    for index in range(len(grouplist)):
                        s += input[index + 2]
                    pval[val].append(round(1-hypergeom.cdf(k-1,M,s,N), 2))
        group += 1
    return pval

def mul_test_corr(list_of_pval, corr_meth, pval):
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
    tests = multitest.multipletests(pvals = list_of_pval, alpha = pval, method = corr_meth)
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
    sig_dif = []
    for TF in TF_list:
        if TF == True:
            sig_dif.append(df_substr.iloc[index]['Substructure'])
        index += 1
    sig_dif.sort(key=len, reverse=True)
    return sig_dif

def create_groups_substr(dif_sig):
    new_dic = {}
    for val in dif_sig:
        placed = 'no'
        for keyval in new_dic.keys():
            if Chem.MolFromSmiles(keyval, sanitize=False).HasSubstructMatch(Chem.MolFromSmiles(val, sanitize=False)):
                new_dic[keyval].append(val)
                placed = 'yes'
        if placed == 'no':
            new_dic[val] = []
    return new_dic
