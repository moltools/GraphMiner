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
#                 print(value)
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

### PREPARATION OF SMILES ###

## Selection on size
from GraphMiner import select_on_size

for group in grouplist:
    list_of_smiles = dict_of_data[group]
    for mol_smile in list_of_smiles:
        select_on_size(mol_smile)

### GRAPH MINING ###

## Breadth First Search, self parsing

from GraphMiner import list_of_nodes, dict_of_nodes, breadth_first_search

#test_smiles = 'NNC(O)=NN'
mol_list = ['C', 'N', 'O', 'Br', 'S']

# nodelist = list_of_nodes(test_smiles, mol_list)
# nodedict = dict_of_nodes(test_smiles, mol_list)
# print(breadth_first_search(nodelist, test_smiles, nodedict))

## Using RDKit to find the neighboring atoms
from GraphMiner import rdkit_parse, list_nodes, breadth_fs, subgraphs_smiles

for group in grouplist:
    list_of_smiles = dict_of_data[group]
    for mol_smile in list_of_smiles:
        list_node = list_nodes(mol_smile, mol_list)
        dictnode = rdkit_parse(mol_smile, list_node)
        subgraphdict = breadth_fs(list_node, dictnode)
        print(subgraphs_smiles(subgraphdict, mol_smile))

        