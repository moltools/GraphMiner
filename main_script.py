import csv


### DATA LOADING ###
from GraphMiner import load_data, determine_groups, create_dict

infile = load_data(1, ';')
grouplist = determine_groups(infile)
dict_of_data = create_dict(grouplist, infile)


### RDKit MINING ###
from GraphMiner import subgraph_miner,sub_to_smiles

### Trying to remove the values that result in errors

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         sub_graphs, m_1, m_ha = subgraph_miner(mol_smile)
#         for list_of_graphs in sub_graphs:
#             for value in list_of_graphs:
#                 for number in value:
#                     if number >= m_ha:
#                         list_of_graphs.remove(value)
#         print(sub_to_smiles(sub_graphs, m_1))

### Basic code that works on simple molecules

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         sub_graphs, m_1, m_ha = subgraph_miner(mol_smile)
#         print(sub_to_smiles(sub_graphs, m_1))

### GRAPH MINING ###

## Breadth First Search, self parsing
from GraphMiner import list_of_nodes, dict_of_nodes, breadth_first_search

# test_smiles = 'NNC(O)=NN'
# mol_list = ['C', 'N', 'O', 'Br', 'S']

# nodelist = list_of_nodes(test_smiles, mol_list)
# nodedict = dict_of_nodes(test_smiles, mol_list)
# print(breadth_first_search(nodelist, test_smiles, nodedict))

## BFS - Using RDKit to find the neighboring atoms without preselection
from GraphMiner import rdkit_parse, breadth_fs, rdkit_smiles

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         dictnode, list_node = rdkit_parse(mol_smile)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         print(subgraphdict)
#         print('Without Preselection')
#         print(subgraphs_smiles(subgraphdict, mol_smile))

## DFS - Using RDKit to find the neighboring atoms without preselection
from GraphMiner import depth_fs

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         dictnode, list_node = rdkit_parse(mol_smile)
#         print(dictnode)
#         print(list_node)
#         dfs_subgraph = depth_fs(list_node, dictnode)
#         print(dfs_subgraph)

### PREPARATION OF SMILES + GRAPH MINING ###

## Selection on size
from GraphMiner import select_on_size

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
        # selected_smile = select_on_size(mol_smile)
        # if selected_smile == None:
        #     continue
#         dictnode, list_node = rdkit_parse(selected_smile)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         # print('With Preselection on Size')
#         smilesdict = rdkit_smiles(subgraphdict, selected_smile)
#         print(' ')
#         print(smilesdict)

## Combination of C-OH, C=O and COOH
from GraphMiner import combine_basic_substructures, return_basic_substructures

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         selected_smile = select_on_size(mol_smile)
#         if selected_smile == None:
#             continue
#         print(' ')
#         print(selected_smile)
#         com_mol_smiles, replaced = combine_basic_substructures(selected_smile)
#         print(com_mol_smiles)
#         dictnode, list_node = rdkit_parse(com_mol_smiles)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         returned_dict = return_basic_substructures(replaced, subgraphdict)
#         print(rdkit_smiles(returned_dict, selected_smile))

# from GraphMiner import combining, returning
# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         selected_smile = select_on_size(mol_smile)
#         if selected_smile == None:
#             continue
#         print(' ')
#         print(selected_smile)
#         com_mol_smiles, new_indeces = combining(selected_smile)
#         print(com_mol_smiles)
#         # dictnode, list_node = rdkit_parse(com_mol_smiles)
#         # subgraphdict = breadth_fs(list_node, dictnode)
#         # print(subgraphdict)
#         repl_smile = returning(com_mol_smiles, new_indeces)
#         print(repl_smile)

## Working with RWMol
from GraphMiner import starting_rwmol

for group in grouplist:
    list_of_smiles = dict_of_data[group]
    for mol_smile in list_of_smiles:
        selected_smile = select_on_size(mol_smile)
        if selected_smile == None:
            continue
        print(' ')
        print(selected_smile)
        print(starting_rwmol(selected_smile))


##Write csv file with frequencies and percentages
from GraphMiner import combine_substr, count_freq, perc_substr

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     all_substr = []
#     total_molecules = 0
#     for mol_smile in list_of_smiles:
#         selected_smile = select_on_size(mol_smile)
#         if selected_smile == None:
#             continue
#         dictnode, list_node = rdkit_parse(selected_smile)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         smilesdict = rdkit_smiles(subgraphdict, selected_smile)
#         unique_str = combine_substr(smilesdict)
#         all_substr += (unique_str)
#         total_molecules += 1
#     print(group)
#     counts = count_freq(all_substr)
#     name = 'Frequency_overview_group' + str(group) +'.csv'
#     percentages = perc_substr(counts, total_molecules)
#     f = open(name,  'w')
#     writer = csv.writer(f)
#     Head_row = ('Substructure', 'Frequency', 'Percentage' + str(group))
#     writer.writerow(Head_row)
#     for perc in percentages:
#         writer.writerow(perc)
#     f.close()


## Write second type CSV file
# from GraphMiner import list_maker
#
# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     all_substr = []
#     dict_substr = {}
#     total_molecules = 0
#     for mol_smile in list_of_smiles:
#         # print('tot: ' + str(total_molecules))
#         selected_smile = select_on_size(mol_smile)
#         if selected_smile == None:
#             continue
#         dictnode, list_node = rdkit_parse(selected_smile)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         smilesdict = rdkit_smiles(subgraphdict, selected_smile)
#         unique_str = combine_substr(smilesdict)
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
