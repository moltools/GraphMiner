#!/usr/bin/env python3

#IMPORTS
from statsmodels.stats import multitest
# import scipy.stats
import pandas as pd

def new_dataframes(input_df):
    input_df.pop('Frequency')
    list_substr = input_df.iloc[:, 0]
    return input_df, list_substr

def join_df(list_of_df):
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    return joined_df

def add_in_perc(joined_substr, all_df, listgroups):
    new_list = [0]*len(joined_substr.iloc[:, 0])
    #WERKT NIET
    for index in range(len(listgroups)):
        joined_substr = joined_substr.assign(percentage = new_list)
    # for substr in joined_substr.iloc[:, 0]:
        # print(len(joined_substr.iloc[:, 0]))
        # print(substr)
    print(joined_substr)
    return new_list


def retrieve_pval():
    return

def bonferonni_corr():
    tests = multitest.multipletests(pvals = [0.1, 0.05, 0.003], alpha = 0.05, method = 'bonferroni')
    return tests

def benj_hoch():
    return
