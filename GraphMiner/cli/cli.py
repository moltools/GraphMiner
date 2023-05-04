#!/usr/bin/env python3
import argparse

def cli():
    parser = argparse.ArgumentParser(prog='GraphMiner', \
                                     description='...',\
                                     epilog='...')
    parser.add_argument("-i", "--input", type=str, default=True, help="input file SMILES and groups, should be a csv file")
    parser.add_argument("-mtc", "--MultipleTestCorr", type=str, choices=['bonferroni', 'fdr_bh', 'holm'], \
                        default='fdr_bh', help="Multiple Testing Correction, choice between bonferroni, fdr_bh and holm")
    parser.add_argument("-o", "--output", type=str, help="Name of the output file")
    parser.add_argument("-size", "--moleculesize", type=int, default=40, help='Number of heavy atoms above which the molecules are not included')
    parser.add_argument("-sep", "--separatorCSVfile", type=str, default=',', help='The separator in the input csv file')
    return parser.parse_args()