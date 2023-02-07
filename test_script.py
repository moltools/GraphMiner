# ###DATA LOADING###
# from GraphMiner import load_data, determine_groups, create_dict

# infile = load_data()
# grouplist = determine_groups(infile)
# dict_of_data = create_dict(grouplist, infile)



###RDKit MINING###
from GraphMiner import subgraph_miner,sub_to_smiles

sub_graphs, m_1 = subgraph_miner('CCC1C(C(C(C(=O)C(CC(C(C(C(C(C(=O)O1)C)OC2CC(C(C(O2)C)O)(C)OC)C)OC3C(C(CC(O3)C)N(C)C)O)(C)O)C)C)O)(C)O')
print(sub_to_smiles(sub_graphs, m_1))