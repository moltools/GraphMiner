from GraphMiner import increment
from GraphMiner import load_data, determine_groups, create_dict

#from GraphMiner import 

infile = load_data()
grouplist = determine_groups(infile)
dict_of_data = create_dict(grouplist, infile)
print((dict_of_data))