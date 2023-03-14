#!/usr/bin/env python3

#IMPORTS
from statsmodels.stats import multitest
import scipy.stats
import pandas as pd
import numpy as np

def new_dataframes(input_df, groupnum):
    '''
    Change imported dataframes to new format for statistical analysis

    input:
    input_df - pandas dataframe containing three columns - substructure in
    SMILES (str), frequency[groupnum] as number of occurrences in int and
    OccurrenceList[groupnum] displaying in str whether the substructure
    occurred (1) or not (0) in a molecule at the index of the molecule
    groupnum - int displaying number of group of molecules

    returns:
    input_df - pandas dataframe containing two columns - substructure in
    SMILES (str) and OccurrenceList[groupnum] displaying whether the
    substructure occurred (1) or not (0) in a molecule at the index of the
    molecule
    list_of_substr - list containing all substructures
    '''
    input_df.pop('Frequency'+str(groupnum))
    list_substr = input_df.iloc[:, 0]
    return input_df, list_substr

def join_df(list_of_df):
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
    #FIND A WAY TO NOT HARDCODE THE ZEROS
    joined_df = joined_df.replace(np.nan, '0,0,0,0')
    return joined_df

def retrieve_pval(df_joined):
    '''
    Calculate the p-values for each susbstructure using two-sided t-test

    input:
    df_joined - dataframe containing three columns - substructure in SMILES
    (str) and twice an OccurrenceList[groupnumber] displaying in str
    whether the substructure occurred (1) or not (0) in a molecule at the
    index of the molecule

    returns:
    pval_list - list containing all p-values calculated for each substructure
    in order of the substructure in the dataframe
    '''
    pval_list = []
    for substr in range(len(df_joined)):
        doc1 = df_joined.iloc[substr, 1].split(',')
        int_doc1 = [int(val) for val in doc1]
        doc2 = df_joined.iloc[substr, 2].split(',')
        int_doc2 = [int(val) for val in doc2]
        pval = scipy.stats.ttest_ind(int_doc1, int_doc2, alternative='two-sided')
        pval_list.append(pval[1])
    return pval_list

def bonferonni_corr(list_of_pval):
    '''
    Perform Multiple Testing Correction using Bonferonni

    input:
    list_of_pval - list containing all p-values calculated for each substructure
    in order of the substructure in the dataframe

    returns:
    tests - result of multiple testing correction as array containing [1] True/
    False, [2] corrected p-values, [3] corrected alpha for Sidak method and
    [4] corrected alpha for Bonferonni method
    '''
    tests = multitest.multipletests(pvals = list_of_pval, alpha = 0.05, method = 'bonferroni')
    return tests

def benj_hoch(list_of_pval):
    '''
    Perform Multiple Testing Correction using Benjamini-Hochberg

    input:
    list_of_pval - list containing all p-values calculated for each substructure
    in order of the substructure in the dataframe

    returns:
    tests - result of multiple testing correction as array containing [1] True/
    False, [2] corrected p-values, [3] corrected alpha for Sidak method and
    [4] corrected alpha for Bonferonni method
    '''
    tests = multitest.multipletests(pvals=list_of_pval, alpha=0.05,
                                    method='fdr_bh')
    return tests
