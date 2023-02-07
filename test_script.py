###DATA LOADING###
from GraphMiner import load_data, determine_groups, create_dict

infile = load_data()
grouplist = determine_groups(infile)
dict_of_data = create_dict(grouplist, infile)
print((dict_of_data))

###