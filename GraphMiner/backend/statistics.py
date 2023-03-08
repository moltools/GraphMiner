#!/usr/bin/env python3

#IMPORTS
import statsmodels

def bonferonni_corr():
    statsmodels.stats.multitest.multipletests(0.05)
    return
