#!/usr/bin/env python3
# from __future__ import annotations
from .cli import cli
import argparse
import pandas as pd
from rdkit import Chem
import numpy as np
import csv
import os
from GraphMiner import determine_groups, create_dict, select_on_size, \
    combine_substr, repl_atommap_COO, set_atommapnum, rdkit_parse_atommap, \
    repl_atommap_CO, repl_atommap_C_O, count_freq, \
    repl_atommap_NCO, repl_atommap_NOCO, select_mol, \
    repl_atommap_POOO, breadth_fs2, return_replaced2, rdkit_smiles2, \
    repl_atommap_SOO, repl_atommap_SOOO, timeout, TimeoutError, \
    repl_atommap_POOOO, hypergeometric_test_pval, \
    mul_test_corr, extract_signif_substr, create_groups_substr, \
    mol_to_fingerprint, plot_dendrogram, create_groups_dendrogram, \
    draw_mol_fig,tanimoto_coefficient
args = cli()

def data_loading(args):
    df = pd.read_csv(args.input, sep=args.separatorCSVfile)
    grouplist = determine_groups(df)
    dictofdata = create_dict(grouplist, df)
    return grouplist, dictofdata

@timeout(args.TimeOutTimer)
def mol_substr_bfs(selected_mol, all_substr, dict_substr, total_molecules):
    # print(' ')
    repl = {}
    sel_smile = Chem.MolToSmiles(selected_mol, kekuleSmiles = True)
    sel_mol = Chem.MolFromSmiles(sel_smile)
    # print('START: ' + sel_smile)
    set_atommapnum(sel_mol)
    tot_mol = sel_mol
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
        sel_mol, repl = repl_atommap_COO(sel_mol, repl)
        # print('C(=O)O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)(O)O')) == True:
        sel_mol, repl = repl_atommap_POOOO(sel_mol, repl)
        # print('P(=O)(O)(O)O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('P(=O)(O)O')) == True:
        sel_mol, repl = repl_atommap_POOO(sel_mol, repl)
        # print('P(=O)(O)O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)O')) == True:
        sel_mol, repl = repl_atommap_SOOO(sel_mol, repl)
        # print('S(=O)(=O)O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('S(=O)(=O)')) == True:
        sel_mol, repl = repl_atommap_SOO(sel_mol, repl)
        # print('S(=O)(=O) done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('N(O)C(=O)')) == True:
        sel_mol, repl = repl_atommap_NOCO(sel_mol, repl)
        # print('NOCO done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('NC=O')) == True:
        sel_mol, repl = repl_atommap_NCO(sel_mol, repl)
        # print('NC=O done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
        sel_mol, repl = repl_atommap_CO(sel_mol, repl)
        # print('CO done')
    if sel_mol.HasSubstructMatch(Chem.MolFromSmiles('C=O')) == True:
        sel_mol, repl = repl_atommap_C_O(sel_mol, repl)
        # print('C=O done')
    dictnode, list_node = rdkit_parse_atommap(sel_mol)
    subgraphdict = breadth_fs2(dictnode, list_node)
    returned_dict = return_replaced2(repl, subgraphdict)
    smilesdict, moldict = rdkit_smiles2(returned_dict, tot_mol, tot_mol)
    unique_str = combine_substr(smilesdict)
    all_substr += (unique_str)
    dict_substr[total_molecules] = unique_str
    total_molecules += 1
    return dict_substr, total_molecules, all_substr

def writesubstrfile(list_of_df, grouplist, dict_of_groups, pathway):
    substrfile = pathway + '/substrfile.csv'
    if len(list_of_df) == 1:
        list_of_df[0].to_csv(substrfile,
                             header=['Substructure', 'Frequency'])
    elif len(list_of_df) == 2:
        joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
        colmn = list(joined_df.columns)
        joined_df[colmn[1]] = joined_df[colmn[1]].replace(np.nan, 0)
        joined_df[colmn[2]] = joined_df[colmn[2]].replace(np.nan, 0)
        headers = ['Substructure']
        for groupname in grouplist:
            headers.append('Frequency' + str(groupname))
        joined_df.to_csv(substrfile, header=headers)
    elif len(list_of_df) >= 3:
        joined_df = pd.merge(list_of_df[0], list_of_df[1], how='outer')
        for dfnum in range(2, len(list_of_df)):
            joined_df = pd.merge(joined_df, list_of_df[dfnum], how='outer')
        colmn = list(joined_df.columns)
        headers = ['Substructure']
        for num_df in range(len(list_of_df)):
            joined_df[colmn[num_df + 1]] = joined_df[
                colmn[num_df + 1]].replace(np.nan, 0)
        for groupname in grouplist:
            headers.append('Frequency' + str(groupname))
        joined_df.to_csv(substrfile, header=headers)
    f = open('datafile.csv', 'w')
    writer = csv.writer(f)
    writer.writerow(dict_of_groups.keys())
    writer.writerow(dict_of_groups.values())
    f.close()
    return

def calculatepval(args, list_of_groups, grouplist, pathway):
    substrfile = pathway + '/substrfile.csv'
    substr_df = pd.read_csv(substrfile, sep=',')
    pvaldict = hypergeometric_test_pval(list_of_groups, substr_df, grouplist)
    for key in pvaldict.keys():
        substr_df[key] = pvaldict[key]
    pvalfile = pathway + '/pvaloverview.csv'
    f = open(pvalfile, 'w')
    writer = csv.writer(f)
    for row in substr_df.iterrows():
        writer.writerow(row)
    f.close()
    return substr_df, pvaldict

def mtc_clustering(pvaldict, substr_df, grouplist, args, pathway):
    filepath1 = pathway + '/significantsubstr.csv'
    f = open(filepath1, 'w')
    writer = csv.writer(f)
    p = 0
    tryoutfail = 0
    imgpath = pathway + '/Images'
    os.mkdir(imgpath)
    for pvallist in pvaldict.values():
        groupname = grouplist[p]
        writer.writerow(['New Group'])
        writer.writerow([p])
        p += 1
        TF_benj_list = mul_test_corr(pvallist, 'fdr_bh', args.PValue)
        substr_df['True/False Benj-Hoch'] = TF_benj_list
        list_sigdif = extract_signif_substr(TF_benj_list, substr_df)
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
        dist_m = np.zeros((len(list_sigdif), len(list_sigdif)))
        for i, fp_i in enumerate(fps):
            for j, fp_j in enumerate(fps):
                if j > i:
                    coef = 1 - tanimoto_coefficient(fp_i, fp_j)
                    dist_m[i, j] = coef
                    dist_m[j, i] = coef
        namefile = pathway + '/Images/' + str(
            groupname) + '_dendrogram.png'
        dendrogram = plot_dendrogram(dist_m, smilessubstr, namefile)
        valueslist = create_groups_dendrogram(dendrogram)
        # writer.writerow(valueslist)
        filepaths = pathway + '/Images/' + str(
            groupname)
        os.mkdir(filepaths)
        draw_mol_fig(valueslist, filepaths)
    f.close()
    return

def main():
    args = cli()
    cur_path = os.getcwd()
    new_path = cur_path + '/GraphMinerResults'
    os.mkdir(new_path)
    # print(new_path)
    group_list, dict_of_data = data_loading(args)
    number = 0
    TimeOut = 0
    list_of_groups = {}
    list_of_df = []
    for group in group_list:
        list_of_smiles = dict_of_data[group]
        all_substr = []
        dict_substr = {}
        group_tot = 0
        for mol_smile in list_of_smiles:
            first_select = select_on_size(mol_smile, args.moleculesize)
            selected_mol = select_mol(first_select)
            if selected_mol == None:
                continue
            elif type(selected_mol) == list:
                # print('list found')
                number += 1
                # print(number)
                for mol in selected_mol:
                    try:
                        dict_substr, group_tot, all_substr = mol_substr_bfs(
                            Chem.MolFromSmiles(mol), all_substr, dict_substr,
                            group_tot)
                    except TimeoutError:
                        # print('timeout')
                        TimeOut += 1
                        continue
            else:
                number += 1
                # print(number)
                try:
                    dict_substr, group_tot, all_substr = mol_substr_bfs(
                        selected_mol, all_substr, dict_substr, group_tot)
                except TimeoutError:
                    # print('timeout')
                    TimeOut += 1
                    continue
        list_of_groups[group] = group_tot
        counts = count_freq(all_substr)
        count_dict = {}
        id = 0
        for substr in counts:
            count_dict[id] = [substr, counts[substr]]
            id += 1
        columnnames = ['Substructure', 'Frequency' + str(group)]
        df = pd.DataFrame.from_dict(count_dict, orient='index',
                                    columns=columnnames)
        list_of_df.append(df)
    print('TimedOutMolecules', TimeOut)
    writesubstrfile(list_of_df, group_list, list_of_groups, new_path)
    substr_df, pvaldict = calculatepval(args, list_of_groups, group_list, new_path)
    mtc_clustering(pvaldict, substr_df, group_list, args, new_path)
    exit(0)


if __name__ == "__main__":
    main()