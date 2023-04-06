#!/usr/bin/env python3

import csv
from rdkit import Chem

### DATA LOADING ###
from GraphMiner import load_data, determine_groups, create_dict

infile = load_data(1, ';')
grouplist = determine_groups(infile)
dict_of_data = create_dict(grouplist, infile)

### PREPARATION OF SMILES + GRAPH MINING ###

from GraphMiner import select_on_size, breadth_fs, rdkit_smiles, \
    combine_substr, repl_atommap_COO, set_atommapnum, rdkit_parse_atommap, \
    return_replaced, repl_atommap_CO, repl_atommap_C_O, count_freq, \
    list_maker, repl_atommap_NCO, repl_atommap_NOCO, select_mol, \
    repl_atommap_POOO

number = 0
# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     all_substr = []
#     dict_substr = {}
#     total_molecules = 0
#     for mol_smile in list_of_smiles:
#         number += 1
#         print(number)
#         first_select = select_on_size(mol_smile)
#         selected_mol = select_mol(first_select)
#         if selected_mol == None:
#             continue
#         print(' ')
#         repl = {}
#         sel_smile = Chem.MolToSmiles(selected_mol, kekuleSmiles = True)
#         sel_mol = Chem.MolFromSmiles(sel_smile)
#         print('START: ' + sel_smile)
#         set_atommapnum(sel_mol)
#         if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
#             sel_mol, repl = repl_atommap_COO(sel_mol, repl)
#             print('C(=O)O done')
#         if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)O')) == True:
#             sel_mol, repl = repl_atommap_POOO(sel_mol, repl)
#             print('P(=O)(O)O done')
#         # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('N(O)C(=O)')) == True:
#         #     sel_mol, repl = repl_atommap_NOCO(sel_mol, repl)
#         #     print('NOCO done')
#         # print(Chem.MolToSmiles(sel_mol))
#         # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('NC=O')) == True:
#         #     sel_mol, repl = repl_atommap_NCO(sel_mol, repl)
#         #     print('NC=O done')
#         if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
#             sel_mol, repl = repl_atommap_CO(sel_mol, repl)
#             print('CO done')
#             # print(Chem.MolToSmiles(sel_mol))
#         if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C=O')) == True:
#             sel_mol, repl = repl_atommap_C_O(sel_mol, repl)
#             print('C=O done')
#         # print(Chem.MolToSmiles(sel_mol))
#         # print(repl)
#         # print(Chem.MolToSmiles(sel_mol))
#         dictnode, list_node = rdkit_parse_atommap(sel_mol)
#         print('parsing done')
#         subgraphdict = breadth_fs(list_node, dictnode)
#         print('bfs done')
#         returned_dict = return_replaced(repl, subgraphdict)
#         # print(returned_dict)
#         print('replaced done')
#         smilesdict = rdkit_smiles(returned_dict, sel_smile)
#         print('smiles returned')
#         # print(smilesdict)
#         unique_str = combine_substr(smilesdict)
#         # print(unique_str)
#         print('all okay')
#         all_substr += (unique_str)
#         dict_substr[total_molecules] = unique_str
#         total_molecules += 1
#     print(group)
#     counts = count_freq(all_substr)
#     list_of_rows = list_maker(dict_substr, counts)
#     name = 'com_overview_group' + str(group) +'.csv'
#     f = open(name,  'w')
#     writer = csv.writer(f)
#     Head_row = ('Substructure', 'Frequency' + str(group), 'OccurrenceList' + str(group))
#     writer.writerow(Head_row)
#     for row in list_of_rows:
#         writer.writerow(row)
#     f.close()
#     print('csv file written')


### STATISTICS PART ###

# Load in csv files
from GraphMiner import new_dataframes

df_list = []
sub_list = []
start = 0
for group in grouplist:
    freq_file = load_data(start+2, ',')
    red_file, substrlist = new_dataframes(freq_file, group)
    df_list.append(red_file)
    start += 1

## Calculate p values
from GraphMiner import join_df, retrieve_pval

substr_df = join_df(df_list)
pvalues = retrieve_pval(substr_df)

## Multiple Testing Correction
from GraphMiner import bonferonni_corr, benj_hoch

TF_bonn_list = bonferonni_corr(pvalues)
TF_benj_list = benj_hoch(pvalues)

#Create CSV file
substr_df['True/False Bonferonni'] = TF_bonn_list
substr_df['True/False Benj-Hoch'] = TF_benj_list
substr_df.to_csv('substructuretruefalse.csv', sep=';')

## Analysis of results
from GraphMiner import extract_signif_substr, create_groups_substr

list_sigdif = extract_signif_substr(TF_benj_list, substr_df)
dic_of_substr = create_groups_substr(list_sigdif)
print(dic_of_substr)






# truetotal = 0
# falsetotal = 0
# for TF in TF_benj_list:
#     if TF == True:
#         truetotal += 1
#     elif TF == False:
#         falsetotal +=1
# print(truetotal, falsetotal)

