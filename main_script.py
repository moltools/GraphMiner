#!/usr/bin/env python3

# LOGGER_LEVEL = 'DEBUG'  # prints debug messages to stdout

import csv
import time
start = time.time()
from rdkit import Chem
import pandas as pd
import numpy as np

# from tqdm import tqdm

### DATA LOADING ###
from GraphMiner import load_data, determine_groups, create_dict

infile = load_data(1, ';')
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


@timeout(30)
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
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)(O)O')) == True:
        sel_mol, repl = repl_atommap_POOOO(sel_mol, repl)
        print('P(=O)(O)(O)O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)O')) == True:
        sel_mol, repl = repl_atommap_POOO(sel_mol, repl)
        print('P(=O)(O)O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)O')) == True:
        sel_mol, repl = repl_atommap_SOOO(sel_mol, repl)
        print('S(=O)(=O)O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)')) == True:
        sel_mol, repl = repl_atommap_SOO(sel_mol, repl)
        print('S(=O)(=O) done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('N(O)C(=O)')) == True:
        sel_mol, repl = repl_atommap_NOCO(sel_mol, repl)
        print('NOCO done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('NC=O')) == True:
        sel_mol, repl = repl_atommap_NCO(sel_mol, repl)
        print('NC=O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
        sel_mol, repl = repl_atommap_CO(sel_mol, repl)
        print('CO done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C=O')) == True:
        sel_mol, repl = repl_atommap_C_O(sel_mol, repl)
        print('C=O done')
    dictnode, list_node = rdkit_parse_atommap(sel_mol)
    subgraphdict = breadth_fs2(dictnode, list_node)
    returned_dict = return_replaced2(repl, subgraphdict)
    smilesdict, moldict = rdkit_smiles2(returned_dict, tot_mol, tot_mol)
    unique_str = combine_substr(smilesdict)
    all_substr += (unique_str)
    dict_substr[total_molecules] = unique_str
    total_molecules += 1
    return dict_substr, total_molecules, all_substr

@timeout(30)
def mol_substr_dfs(selected_mol, all_substr, dict_substr, total_molecules):
    print(' ')
    repl = {}
    sel_smile = Chem.MolToSmiles(selected_mol, kekuleSmiles = True)
    sel_mol = Chem.MolFromSmiles(sel_smile)
    print('START: ' + sel_smile)
    set_atommapnum(sel_mol)
    tot_mol = sel_mol
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
    #     sel_mol, repl = repl_atommap_COO(sel_mol, repl)
    #     print('C(=O)O done')
    #     # print(Chem.MolToSmiles(sel_mol))
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)O')) == True:
    #     sel_mol, repl = repl_atommap_POOO(sel_mol, repl)
    #     print('P(=O)(O)O done')
    #     # print(Chem.MolToSmiles(sel_mol))
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)O')) == True:
    #     sel_mol, repl = repl_atommap_SOOO(sel_mol, repl)
    #     print('S(=O)(=O)O done')
    #     # print(Chem.MolToSmiles(sel_mol))
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)')) == True:
    #     sel_mol, repl = repl_atommap_SOO(sel_mol, repl)
    #     print('S(=O)(=O) done')
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('N(O)C(=O)')) == True:
    #     sel_mol, repl = repl_atommap_NOCO(sel_mol, repl)
    #     print('NOCO done')
    #     # print(Chem.MolToSmiles(sel_mol))
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('NC=O')) == True:
    #     sel_mol, repl = repl_atommap_NCO(sel_mol, repl)
    #     print('NC=O done')
    #     # print(Chem.MolToSmiles(sel_mol))
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
    #     sel_mol, repl = repl_atommap_CO(sel_mol, repl)
    #     print('CO done')
    #     # print(Chem.MolToSmiles(sel_mol))
    # if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C=O')) == True:
    #     sel_mol, repl = repl_atommap_C_O(sel_mol, repl)
    #     print('C=O done')
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
    return dict_substr, total_molecules, all_substr

number = 0
group_num = 0
len_dict = {}
list_of_groups = {}
list_of_df = []
TO = 0
TOlist = []
tottlist = []
for group in grouplist:
    list_of_smiles = dict_of_data[group]
    all_substr = []
    dict_substr = {}
    group_tot = 0
    for mol_smile in list_of_smiles:
        first_select = select_on_size(mol_smile, 40)
        selected_mol = select_mol(first_select)
        if selected_mol == None:
            continue
        elif type(selected_mol) == list:
            print('list found')
            number += 1
            print(number)
            for mol in selected_mol:
                try:
                    dict_substr, group_tot, all_substr = mol_substr_bfs(Chem.MolFromSmiles(mol), all_substr, dict_substr, group_tot)
                    tottlist.append(Chem.MolFromSmiles(mol).GetNumHeavyAtoms())
                except TimeoutError:
                    print('timeout')
                    TO += 1
                    TOlist.append(Chem.MolFromSmiles(mol).GetNumHeavyAtoms())
                    continue
        else:
            number += 1
            print(number)
            try:
                dict_substr, group_tot, all_substr = mol_substr_bfs(selected_mol, all_substr, dict_substr, group_tot)
                tottlist.append(selected_mol.GetNumHeavyAtoms())
            except TimeoutError:
                print('timeout')
                TO += 1
                TOlist.append(selected_mol.GetNumHeavyAtoms())
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
print('TO', TO)
print(TOlist)
print(tottlist)
if len(list_of_df) == 1:
    list_of_df[0].to_csv('substrfile.csv', header = ['Substructure', 'Frequency'])
elif len(list_of_df) == 2:
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    colmn = list(joined_df.columns)
    joined_df[colmn[1]] = joined_df[colmn[1]].replace(np.nan, 0)
    joined_df[colmn[2]] = joined_df[colmn[2]].replace(np.nan, 0)
    headers = ['Substructure']
    for groupname in grouplist:
        headers.append('Frequency' + str(groupname))
    print(headers)
    joined_df.to_csv('substrfile.csv', header = headers)
elif len(list_of_df) >=3:
    joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
    for dfnum in range(2, len(list_of_df)):
        joined_df = pd.merge(joined_df, list_of_df[dfnum], how='outer')
    print(joined_df)
    colmn = list(joined_df.columns)
    headers = ['Substructure']
    for num_df in range(len(list_of_df)):
        joined_df[colmn[num_df + 1]] = joined_df[colmn[num_df + 1]].replace(np.nan, 0)
    for groupname in grouplist:
        headers.append('Frequency' + str(groupname))
    print(headers)
    joined_df.to_csv('substrfile.csv',header=headers)

f = open('datafile.csv', 'w')
writer = csv.writer(f)
writer.writerow(grouplist)
writer.writerow(list_of_groups.keys())
writer.writerow(list_of_groups.values())
f.close()

### STATISTICS PART ###

# Load in csv files
substr_df = load_data(2, ',')

## Calculate p values
from GraphMiner import hypergeometric_test_pval, \
    mul_test_corr, extract_signif_substr, create_groups_substr, \
    mol_to_fingerprint, plot_dendrogram, create_groups_dendrogram, \
    draw_mol_fig,tanimoto_coefficient

pvaldict = hypergeometric_test_pval(list_of_groups, substr_df, grouplist)
for key in pvaldict.keys():
    substr_df[key] = pvaldict[key]

f = open('pvaloverview.csv',  'w')
writer = csv.writer(f)
for row in substr_df.iterrows():
    writer.writerow(row)
f.close()

f = open('significantsubstr.csv', 'w')
writer = csv.writer(f)
p = 0
tryoutfail = 0
groupnum = 0
for pvallist in pvaldict.values():
    writer.writerow(['New Group'])
    writer.writerow([p])
    p +=1
    TF_benj_list = mul_test_corr(pvallist, 'fdr_bh', 0.05)
    substr_df['True/False Benj-Hoch'] = TF_benj_list
    list_sigdif = extract_signif_substr(TF_benj_list, substr_df)
    # writer.writerow(list_sigdif)
    dic_of_substr = create_groups_substr(list_sigdif)
    writer.writerow(dic_of_substr.keys())
    if len(list_sigdif) < 3:
        continue
    smilessubstr = []
    for smiless in list_sigdif:
        smiless = smiless.replace('c', 'C')
        smiless = smiless.replace('n', 'N')
        smiless = smiless.replace('o', 'O')
        smiless = smiless.replace('s', 'S')
        smilessubstr.append(smiless)
    molsubstr = [Chem.MolFromSmiles(smiles) for smiles in smilessubstr]
    fps = []
    for mol in molsubstr:
        try:
            fpsmol = mol_to_fingerprint(mol)
            fps.append(fpsmol)
        except:
            tryoutfail += 1
            # logger.debug(Exception)
            # logger.debug(type(Exception).__name__)
            continue
            # fps = [mol_to_fingerprint(mol) for mol in molsubstr]
    dist_m = np.zeros((len(list_sigdif), len(list_sigdif)))
    for i, fp_i in enumerate(fps):
        for j, fp_j in enumerate(fps):
            if j > i:
                coef = 1 - tanimoto_coefficient(fp_i, fp_j)
                dist_m[i, j] = coef
                dist_m[j, i] = coef
    print(dist_m)
    namefile = '/home/duive014/MOLTOOLS/GraphMiner/Images/' + str(groupnum) + '_dendrogram.png'
    dendrogram = plot_dendrogram(dist_m, smilessubstr, namefile)
    print(groupnum)
    valueslist = create_groups_dendrogram(dendrogram)
    # writer.writerow(valueslist)
    filepaths = '/home/duive014/MOLTOOLS/GraphMiner/Images' + '/' + str(groupnum)
    draw_mol_fig(valueslist, filepaths)
    groupnum += 1
f.close()
end = time.time()
print('time', end-start)
print('TO', TO)
print('number', number)
print('tryoutfail', tryoutfail)


# if __name__ == "__main__":
#     logger = Logger('main_script.py')

#     stream_handler = StreamHandler()
#     stream_handler.setLevel(LOGGER_LEVEL)
#     logger.addHandler(stream_handler)

#     file_handler = FileHandler('main_script.log')  # this is the file where all debug messages will be written to
#     file_handler.setLevel("DEBUG")  # always log debug messages to file
#     logger.addHandler(file_handler)

#     main(logger)

