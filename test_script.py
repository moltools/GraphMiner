### DATA LOADING ###
# from GraphMiner import load_data, determine_groups, create_dict

# infile = load_data()
# grouplist = determine_groups(infile)
# dict_of_data = create_dict(grouplist, infile)


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

### GRAPH MINING ###
from GraphMiner import list_of_nodes, dict_of_nodes, breadth_first_search

mol_list = ['C', 'N', 'O', 'Br']
nodelist = list_of_nodes('CC(=O)C(=N)Br', mol_list)
print(breadth_first_search(nodelist, 'CC(=O)C(=N)Br'))