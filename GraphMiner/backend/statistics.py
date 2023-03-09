#!/usr/bin/env python3

#IMPORTS
from statsmodels.stats import multitest
# import scipy.stats

def new_dataframes(input_df):
    input_df.pop('Frequency')
    return input_df

def retrieve_pval():
    return

def bonferonni_corr():
    tests = multitest.multipletests(pvals = [0.1, 0.05, 0.003], alpha = 0.05, method = 'bonferroni')
    return tests

def benj_hoch():
    return
