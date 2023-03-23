#!/usr/bin/env python3

import csv
from rdkit import Chem

### DATA LOADING ###
from GraphMiner import load_data, determine_groups, create_dict

infile = load_data(1, ',')
grouplist = determine_groups(infile)
dict_of_data = create_dict(grouplist, infile)

### PREPARATION OF SMILES + GRAPH MINING ###

## Working with RWMol
from GraphMiner import select_on_size, breadth_fs, rdkit_smiles, combine_substr, \
    repl_atommap_COO, set_atommapnum, rdkit_parse_atommap, return_replaced, repl_atommap_CO, \
    repl_atommap_C_O

for group in grouplist:
    list_of_smiles = dict_of_data[group]
    for mol_smile in list_of_smiles:
        selected_mol = select_on_size(mol_smile)
        if selected_mol == None:
            continue
        print(' ')
        repl = {}
        sel_smile = Chem.MolToSmiles(selected_mol)
        sel_mol = Chem.MolFromSmiles(sel_smile)
        print('START2: ' + sel_smile)
        set_atommapnum(sel_mol)
        if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
            sel_mol, repl = repl_atommap_COO(sel_mol, repl)
            print('C(=O)O done')
        # print(Chem.MolToSmiles(sel_mol))
        if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
            sel_mol, repl = repl_atommap_CO(sel_mol, repl)
            print('CO done')
        # print(Chem.MolToSmiles(sel_mol))
        if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C=O')) == True:
            sel_mol, repl = repl_atommap_C_O(sel_mol, repl)
            print('C=O done')
        # print(repl)
        dictnode, list_node = rdkit_parse_atommap(sel_mol)
        subgraphdict = breadth_fs(list_node, dictnode)
        returned_dict = return_replaced(repl, subgraphdict)
        # print(returned_dict)
        smilesdict = rdkit_smiles(returned_dict, sel_smile)
        # print(smilesdict)
        unique_str = combine_substr(smilesdict)
        # print(unique_str)
        print('all okay')

# from GraphMiner import combining, returning
# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         selected_smile = select_on_size(mol_smile)
#         if selected_smile == None:
#             continue
#         print(' ')
#         com_mol_smiles, new_indeces = combining(Chem.MolToSmiles(selected_smile))
#         print(com_mol_smiles)
#         dictnode, list_node = rdkit_parse(com_mol_smiles)
#         # subgraphdict = breadth_fs(list_node, dictnode)
#         # print(subgraphdict)
#         repl_smile = returning(com_mol_smiles, new_indeces)
#         print(repl_smile)


## Write second type CSV file
from GraphMiner import list_maker

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     all_substr = []
#     dict_substr = {}
#     total_molecules = 0
#     for mol_smile in list_of_smiles:
#         # print('tot: ' + str(total_molecules))
#         selected_mol = select_on_size(mol_smile)
#         if selected_mol == None:
#             continue
#         selected_smile = Chem.MolToSmiles(selected_mol)
#         dictnode, list_node = rdkit_parse(selected_smile)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         smilesdict = rdkit_smiles(subgraphdict, selected_smile)
#         unique_str = combine_substr(smilesdict)
#         print(unique_str)
#         all_substr += (unique_str)
#         dict_substr[total_molecules] = unique_str
#         total_molecules += 1
#     print(group)
#     counts = count_freq(all_substr)
#     list_of_rows = list_maker(dict_substr, counts)
#     name = 'New_overview_group' + str(group) +'.csv'
#     f = open(name,  'w')
#     writer = csv.writer(f)
#     Head_row = ('Substructure', 'Frequency' + str(group), 'OccurrenceList' + str(group))
#     writer.writerow(Head_row)
#     for row in list_of_rows:
#         writer.writerow(row)
#     f.close()

### STATISTICS PART ###

## Load in csv files
from GraphMiner import new_dataframes

# df_list = []
# sub_list = []
# for group in grouplist:
#     freq_file = load_data(group+2, ',')
#     red_file, substrlist = new_dataframes(freq_file, group)
#     df_list.append(red_file)

## Calculate p values
from GraphMiner import join_df, retrieve_pval

# substr_df = join_df(df_list)
# pvalues = retrieve_pval(substr_df)

## Multiple Testing Correction
from GraphMiner import bonferonni_corr, benj_hoch
#
# mtc_bonn = bonferonni_corr(pvalues)
# TF_bonn_list = mtc_bonn[0].tolist()
# mtc_benj = benj_hoch(pvalues)
# TF_benj_list = mtc_benj[0].tolist()
# substr_df['True/False Bonferonni'] = TF_bonn_list
# substr_df['True/False Benj-Hoch'] = TF_benj_list
# substr_df.to_csv('substructuretruefalse.csv', sep=';')
#
# truetotal = 0
# falsetotal = 0
# for TF in TF_benj_list:
#     if TF ==True:
#         truetotal += 1
#     elif TF == False:
#         falsetotal +=1
# print(truetotal, falsetotal)
