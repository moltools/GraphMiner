#!/usr/bin/env python3
from .cli import cli
import argparse
from .backend import determine_groups, create_dict, select_on_size, select_mol
import pandas as pd
from rdkit import Chem

def main():
    args = cli()
    # print('Hi Giovi')

    ### DATA LOADING ###
    df = pd.read_csv(args.input, args.separatorCSVfile)
    grouplist = determine_groups(df)
    dict_of_data = create_dict(grouplist, df)

    ### PREPARATION OF SMILES + GRAPH MINING ###
    for group in grouplist:
        list_of_smiles = dict_of_data[group]
        all_substr = []
        dict_substr = {}
        total_molecules = 0
        for mol_smile in list_of_smiles:
            first_select = select_on_size(mol_smile, args.moleculesize)
            selected_mol = select_mol(first_select)
            if selected_mol == None:
                continue
            print(Chem.MolToSmiles(selected_mol))
    exit(0)


if __name__ == "__main__":
    main()