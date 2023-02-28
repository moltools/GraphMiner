import csv

### DATA LOADING ###
from GraphMiner import load_data, determine_groups, create_dict

infile = load_data()
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

test_smiles = 'NNC(O)=NN'
mol_list = ['C', 'N', 'O', 'Br', 'S']

# nodelist = list_of_nodes(test_smiles, mol_list)
# nodedict = dict_of_nodes(test_smiles, mol_list)
# print(breadth_first_search(nodelist, test_smiles, nodedict))

## BFS - Using RDKit to find the neighboring atoms without preselection
from GraphMiner import rdkit_parse, breadth_fs, rdkit_smiles

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         dictnode, list_node = rdkit_parse(selected_smile)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         print(subgraphdict)
#         print('Without Preselection')
#         print(subgraphs_smiles(subgraphdict, mol_smile))

## DFS - Using RDKit to find the neighboring atoms without preselection
from GraphMiner import depth_fs

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         dictnode, list_node = rdkit_parse(selected_smile)
#         dfs_subgraph = depth_fs(list_node, dictnode)
#         print(dfs_subgraph)

### PREPARATION OF SMILES + GRAPH MINING ###

## Selection on size
from GraphMiner import select_on_size

# for group in grouplist:
#     list_of_smiles = dict_of_data[group]
#     for mol_smile in list_of_smiles:
#         selected_smile = select_on_size(mol_smile)
#         if selected_smile == None:
#             continue
#         dictnode, list_node = rdkit_parse(selected_smile)
#         subgraphdict = breadth_fs(list_node, dictnode)
#         print('With Preselection on Size')
#         smilesdict = rdkit_smiles(subgraphdict, selected_smile)
#         print(smilesdict)

## Combination of C-OH, C=O and COOH
from GraphMiner import combine_basic_substructures, return_basic_substructures

for group in grouplist:
    list_of_smiles = dict_of_data[group]
    for mol_smile in list_of_smiles:
        selected_smile = select_on_size(mol_smile)
        if selected_smile == None:
            continue
        com_mol_smiles, replaced = combine_basic_substructures(mol_smile)
        dictnode, list_node = rdkit_parse(com_mol_smiles)
        subgraphdict = breadth_fs(list_node, dictnode)
        print(subgraphdict)
        print(replaced)
        print('With Preselection on Size and Combination')

        returned_dict = return_basic_substructures(selected_smile, replaced, subgraphdict)
        print(rdkit_smiles(returned_dict, selected_smile))

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
#     counts = count_freq(all_substr)
#     name = 'Frequency_overview_group' + str(group) +'.csv'
#     percentages = perc_substr(counts, total_molecules)
#     f = open(name,  'w')
#     writer = csv.writer(f)
#     Head_row = ('Substructure', 'Frequency', 'Percentage')
#     writer.writerow(Head_row)
#     for perc in percentages:
#         writer.writerow(perc)
#     f.close()
    