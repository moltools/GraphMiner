#!/usr/bin/env python3

#IMPORTS
from statsmodels.stats import multitest
import scipy.stats
import pandas as pd
import numpy as np

def new_dataframes(input_df, groupnum):
    input_df.pop('Frequency'+str(groupnum))
    list_substr = input_df.iloc[:, 0]
    return input_df, list_substr

def join_df(list_of_df):
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    #FIND A WAY TO NOT HARDCODE THE ZEROS
    joined_df = joined_df.replace(np.nan, '0,0,0,0')
    return joined_df

def retrieve_pval(df_joined):
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
    tests = multitest.multipletests(pvals = list_of_pval, alpha = 0.05, method = 'bonferroni')
    return tests

def benj_hoch():
    return
