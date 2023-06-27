#!/usr/bin/env python3
import argparse

def cli():
    parser = argparse.ArgumentParser(prog='GraphMiner', \
                                     description='GraphMiner is a tool which performs subgraph mining and enrichment testing on groups of molecules',\
                                     epilog='All output of GraphMiner can be found in one folder, by default named GraphMinerResults.')

    parser.add_argument("-i", "--input", type=str, default=True, \
                        help="input file SMILES and groups, should be a csv file", \
                        required=True)

    parser.add_argument("-mtc", "--MultipleTestCorr", type=str, \
                        choices=['bonferroni', 'fdr_bh', 'holm'], \
                        default='fdr_bh', \
                        help="Multiple Testing Correction, choice between bonferroni, fdr_bh and holm")

    parser.add_argument("-o", "--outputfolder", type=str, \
                        default='GraphMinerResults', \
                        help="Name of the output folder containing all output files")

    parser.add_argument("-size", "--moleculesize", type=int, default=60, \
                        help='Number of heavy atoms above which the molecules are not included')

    parser.add_argument("-sep", "--separatorCSVfile", type=str, default=',', \
                        help='The separator in the input csv file')

    parser.add_argument('-time', '--TimeOutTimer', type=int, default=30, \
                        help='Number of seconds after which the substructure search for a single molecule stops')

    parser.add_argument('-pval', '--PValue', type=float, default=0.05, \
                        help='p-value used for hypergeometric test')

    parser.add_argument('-dend', '--CutOffDendrogram', type=float, default = 1.5, \
                        help='Cut Off value in the Dendrogram to create the groups of under/over enriched molecules')

    return parser.parse_args()