#!/usr/bin/env python3
# from __future__ import annotations
from .cli import cli
import argparse
from .backend import determine_groups, create_dict, select_on_size, \
    select_mol, timeout, set_atommapnum, repl_atommap_COO, repl_atommap_SOO, \
    repl_atommap_SOOO, repl_atommap_POOO, repl_atommap_C_O, repl_atommap_CO, \
    repl_atommap_NCO, repl_atommap_NOCO, rdkit_parse_atommap, breadth_fs2, \
    return_replaced2, rdkit_smiles2, combine_substr
import pandas as pd
from rdkit import Chem

def data_loading(args):
    df = pd.read_csv(args.input, args.separatorCSVfile)
    grouplist = determine_groups(df)
    dictofdata = create_dict(grouplist, df)
    return grouplist, dictofdata

@timeout(5)
def mol_substr_bfs(selected_mol, all_substr, dict_substr, total_molecules):
    print(' ')
    repl = {}
    sel_smile = Chem.MolToSmiles(selected_mol, kekuleSmiles = True)
    sel_mol = Chem.MolFromSmiles(sel_smile)
    print('START: ' + sel_smile)
    set_atommapnum(sel_mol)
    tot_mol = sel_mol
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
        sel_mol, repl = repl_atommap_COO(sel_mol, repl)
        print('C(=O)O done')
        # print(Chem.MolToSmiles(sel_mol))
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)O')) == True:
        sel_mol, repl = repl_atommap_POOO(sel_mol, repl)
        print('P(=O)(O)O done')
        # print(Chem.MolToSmiles(sel_mol))
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)O')) == True:
        sel_mol, repl = repl_atommap_SOOO(sel_mol, repl)
        print('S(=O)(=O)O done')
        # print(Chem.MolToSmiles(sel_mol))
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)')) == True:
        sel_mol, repl = repl_atommap_SOO(sel_mol, repl)
        print('S(=O)(=O) done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('N(O)C(=O)')) == True:
        sel_mol, repl = repl_atommap_NOCO(sel_mol, repl)
        print('NOCO done')
        # print(Chem.MolToSmiles(sel_mol))
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('NC=O')) == True:
        sel_mol, repl = repl_atommap_NCO(sel_mol, repl)
        print('NC=O done')
        # print(Chem.MolToSmiles(sel_mol))
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
        sel_mol, repl = repl_atommap_CO(sel_mol, repl)
        print('CO done')
        # print(Chem.MolToSmiles(sel_mol))
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C=O')) == True:
        sel_mol, repl = repl_atommap_C_O(sel_mol, repl)
        print('C=O done')
        # print(Chem.MolToSmiles(sel_mol))
    dictnode, list_node = rdkit_parse_atommap(sel_mol)
    subgraphdict = breadth_fs2(dictnode, list_node)
    returned_dict = return_replaced2(repl, subgraphdict)
    smilesdict = rdkit_smiles2(returned_dict, tot_mol)
    unique_str = combine_substr(smilesdict)
    all_substr += (unique_str)
    dict_substr[total_molecules] = unique_str
    total_molecules += 1
    return dict_substr, total_molecules

def main():
    args = cli()
    group_list, dict_of_data = data_loading(args)
    if args.SelectionMethod == 'MoleculeSize':
        print('molsize')
        ### PREPARATION OF SMILES + GRAPH MINING ###
        number = 0
        for group in group_list:
            list_of_smiles = dict_of_data[group]
            all_substr = []
            dict_substr = {}
            total_molecules = 0
            Time_Out_Error = 0
            TOlist = []
            for mol_smile in list_of_smiles:
                first_select = select_on_size(mol_smile, args.moleculesize)
                selected_mol = select_mol(first_select)
                if selected_mol == None:
                    continue
                print(Chem.MolToSmiles(selected_mol))
                try:
                    mol_substr_bfs(selected_mol, all_substr, dict_substr,
                               total_molecules)
                except:
                    Time_Out_Error += 1
                    print('Time out', Time_Out_Error)
                    TOlist.append(Chem.MolToSmiles(selected_mol))
                    continue
                finally:
                    number += 1
                    print('number', number)
    elif args.SelectionMethod == 'TimeOutTimer':
        print('TimeOutTimer')
        number = 0
        for group in group_list:
            list_of_smiles = dict_of_data[group]
            all_substr = []
            dict_substr = {}
            total_molecules = 0
            Time_Out_Error = 0
            TOlist = []
            for mol_smile in list_of_smiles:
                selected_mol = select_mol(mol_smile)
                if selected_mol == None:
                    continue
                print(Chem.MolToSmiles(selected_mol))
                try:
                    mol_substr_bfs(selected_mol, all_substr, dict_substr,
                               total_molecules)
                except:
                    Time_Out_Error += 1
                    print('Time out', Time_Out_Error)
                    TOlist.append(Chem.MolToSmiles(selected_mol))
                    continue
                finally:
                    number += 1
                    print('number', number)
    exit(0)


if __name__ == "__main__":
    main()