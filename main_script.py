#!/usr/bin/env python3

import csv
import time
start = time.time()
from rdkit import Chem
import pandas as pd
import numpy as np


### DATA LOADING ###
from GraphMiner import load_data, determine_groups, create_dict

infile = load_data(1, ',')
grouplist = determine_groups(infile)
dict_of_data = create_dict(grouplist, infile)

### PREPARATION OF SMILES + GRAPH MINING ###

from GraphMiner import select_on_size, \
    combine_substr, repl_atommap_COO, set_atommapnum, rdkit_parse_atommap, \
    repl_atommap_CO, repl_atommap_C_O, count_freq, \
    repl_atommap_NCO, repl_atommap_NOCO, select_mol, \
    repl_atommap_POOO, breadth_fs2, depth_fs, return_replaced2, \
    rdkit_smiles2, repl_atommap_SOO, repl_atommap_SOOO, timeout, return_replaced3, \
    rdkit_smiles3, combine_substr2, TimeoutError, repl_atommap_POOOO


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
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)(O)O')) == True:
        sel_mol, repl = repl_atommap_POOOO(sel_mol, repl)
        print('P(=O)(O)(O)O done')
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
    return dict_substr, total_molecules, all_substr

@timeout(120)
def mol_substr_dfs(selected_mol, all_substr, dict_substr, total_molecules):
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
    subgraphdict = depth_fs(sel_mol, dictnode)
    returned_dict = return_replaced3(repl, subgraphdict)
    # print(returned_dict)
    smilesdict = rdkit_smiles3(returned_dict, tot_mol)
    unique_str = combine_substr2(smilesdict)
    all_substr += (unique_str)
    dict_substr[total_molecules] = unique_str
    total_molecules += 1
    return dict_substr, total_molecules

number = 0
group_num = 0
len_dict = {}
list_of_groups = {}
list_of_df = []
for group in grouplist:
    list_of_smiles = dict_of_data[group]
    all_substr = []
    dict_substr = {}
    group_tot = 0
    for mol_smile in list_of_smiles:
        first_select = select_on_size(mol_smile, 10000000)
        selected_mol = select_mol(first_select)
        if selected_mol == None:
            continue
        elif type(selected_mol) == list:
            print('list found')
            number += 1
            for mol in selected_mol:
                try:
                    dict_substr, group_tot, all_substr = mol_substr_bfs(Chem.MolFromSmiles(mol), all_substr, dict_substr, group_tot)
                except TimeoutError:
                    print('timeout')
                    continue
        else:
            number += 1
            try:
                dict_substr, group_tot, all_substr = mol_substr_bfs(selected_mol, all_substr, dict_substr, group_tot)
            except TimeoutError:
                print('timeout')
                continue
    list_of_groups[group] = group_tot
    counts = count_freq(all_substr)
    count_dict = {}
    id=0
    for substr in counts:
        count_dict[id] = [substr, counts[substr]]
        id+=1
    columnnames = ['Substructure', 'Frequency' + str(group_num)]
    df = pd.DataFrame.from_dict(count_dict, orient='index', columns=columnnames)
    print('df')
    print(df)
    list_of_df.append(df)
    group_num += 1
    print(list_of_df)

if len(list_of_df) == 1:
    list_of_df[0].to_csv('substrfile.csv', header = ['Substructure', 'Frequency'])
elif len(list_of_df) == 2:
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    colmn = list(joined_df.columns)
    joined_df[colmn[1]] = joined_df[colmn[1]].replace(np.nan, 0)
    joined_df[colmn[2]] = joined_df[colmn[2]].replace(np.nan, 0)
    joined_df.to_csv('substrfile.csv', header = ['Substructure', 'Frequency0', 'Frequency1'])
elif len(list_of_df) >=3:
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    for dfnum in range(2, len(list_of_df)):
        joined_df = pd.merge(joined_df, list_of_df[dfnum], how='outer')
    print(joined_df)
    colmn = list(joined_df.columns)
    headers = ['Substructure']
    for num_df in range(len(list_of_df)):
        joined_df[colmn[num_df + 1]] = joined_df[colmn[num_df + 1]].replace(np.nan, 0)
        headers.append('Frequency' + str(num_df))
    print(headers)
    joined_df.to_csv('substrfile.csv',header=headers)


### STATISTICS PART ###

# Load in csv files
substr_df = load_data(2, ',')

## Calculate p values
from GraphMiner import hypergeometric_test_pval, \
    mul_test_corr, extract_signif_substr, create_groups_substr

pvaldict = hypergeometric_test_pval(list_of_groups, substr_df, grouplist)
for key in pvaldict.keys():
    substr_df[key] = pvaldict[key]

f = open('pvaloverview.csv',  'w')
writer = csv.writer(f)
for row in substr_df.iterrows():
    writer.writerow(row)
f.close()

for pvallist in pvaldict.values():
    TF_benj_list = mul_test_corr(pvallist, 'fdr_bh')
    substr_df['True/False Benj-Hoch'] = TF_benj_list
    list_sigdif = extract_signif_substr(TF_benj_list, substr_df)
    print(list_sigdif)
    dic_of_substr = create_groups_substr(list_sigdif)
    print(dic_of_substr)

end = time.time()
print('time', end-start)
