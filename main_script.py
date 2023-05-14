#!/usr/bin/env python3

import csv
import time
start = time.time()
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
    repl_atommap_POOO, breadth_fs2, depth_fs, return_replaced2, \
    rdkit_smiles2, repl_atommap_SOO, repl_atommap_SOOO, repl_atommap_benzene, \
    repl_atommap_cyclohex, remove_atom_charges, timeout, return_replaced3, \
    rdkit_smiles3, combine_substr2, list_maker2


@timeout(120)
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
Time_Out_Error = 0
TOlist = []
len_dict = {}
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
        number += 1
        try:
            dict_substr, group_tot = mol_substr_bfs(selected_mol, all_substr, dict_substr, group_tot)
        except:
            # Time_Out_Error +=1
            # print('Time out', Time_Out_Error)
            # TOlist.append(Chem.MolToSmiles(selected_mol))
            continue

    print(group)
    counts = count_freq(all_substr)
    list_of_rows = list_maker2(counts)
    name = 'new_type_csv' + str(group) +'.csv'
    f = open(name,  'w')
    writer = csv.writer(f)
    Head_row = ('Substructure', 'Frequency' + str(group))
    writer.writerow(Head_row)
    for row in list_of_rows:
        writer.writerow(row)
    f.close()
    print('csv file written')


# print(TOlist)
# print(Time_Out_Error)


### STATISTICS PART ###

# Load in csv files
from GraphMiner import new_dataframes

df_list = []
sub_list = []
start = 0
for group in grouplist:
    freq_file = load_data(start+2, ',')
    df_list.append(freq_file)
    start += 1
    print(group)

## Calculate p values
from GraphMiner import join_df2, retrieve_pval

substr_df = join_df2(df_list)
print('dataframe made')
print(substr_df)
##MOET NOG DE TOTALE AANTAL MOLECULEN IN EEN GROEP EN IN TOTAAL STOREN ERGENS ##HAHAH OEEPSIEFLOEPSIE
###MISSCHIEN APARTE CSV FILE HIERVOOR??
# pvalues = retrieve_pval(substr_df)
# print('pvalues calculated')
#
# ## Multiple Testing Correction
# from GraphMiner import mul_test_corr
#
# # TF_bonn_list = mul_test_corr(pvalues, 'bonferroni')
# TF_benj_list = mul_test_corr(pvalues, 'fdr_bh')
# print('mtc done')
#
# #Create CSV file
# # substr_df['True/False Bonferonni'] = TF_bonn_list
# substr_df['True/False Benj-Hoch'] = TF_benj_list
# # substr_df.to_csv('substructuretruefalse.csv', sep=';')
# print('TFlist')
#
# ## Analysis of results
# from GraphMiner import extract_signif_substr, create_groups_substr
#
# list_sigdif, dict_sigdif = extract_signif_substr(TF_benj_list, substr_df)
# print(dict_sigdif)
# dic_of_substr = create_groups_substr(list_sigdif)
# print(dic_of_substr)






# truetotal = 0
# falsetotal = 0
# for TF in TF_benj_list:
#     if TF == True:
#         truetotal += 1
#     elif TF == False:
#         falsetotal +=1
# print(truetotal, falsetotal)

end = time.time()
print('time', end-start)