#!/usr/bin/env python3

#IMPORTS
# from statsmodels.stats import multitest
import scipy.stats
import pandas as pd
import numpy as np

def new_dataframes(input_df):
    input_df.pop('Frequency')
    list_substr = input_df.iloc[:, 0]
    return input_df, list_substr

def join_df(list_of_df):
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    joined_df = joined_df.fillna(0.0)
    return joined_df

def retrieve_pval(df_joined):
    scipy.stats.ttest_ind(df_joined.iloc[:, 1], df_joined.iloc[:, 2], alternative='two-sided')
    return

def bonferonni_corr():
    # tests = multitest.multipletests(pvals = [0.1, 0.05, 0.003], alpha = 0.05, method = 'bonferroni')
    return 

def benj_hoch():
    return
