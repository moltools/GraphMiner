#!/usr/bin/env python3

def breadth_fs(nodelist:list, neigh_dict:dict):
    '''
    Find all subgraphs that are present
    
    input:
    nodelist - list containing the indeces (int) in the string containing an atom
    neigh_dicht - dictionary containing as key(int) the index of the start atom
    and as value (list of int) the indices of the neighbouring atoms

    returns:
    subgraphs - dictionary containg as key (int) the subgraph length and as value
    (list of str) with the subgraphs as strings of indeces divided by '-'
    '''
    len_subgraph = 1
    subgraphs = {}
    subgraphs[len_subgraph] = []
    #Generate the first list of keys
    for key in neigh_dict:
        subgraphs[len_subgraph].append(str(key))
    #Add in all unqiue regions for each length of subgraphs
    for index in range(len(nodelist)):
        len_subgraph += 1
        subgraphs[len_subgraph] = []
        sorted_graphs = []
        for values in subgraphs[len_subgraph-1]:
            for added in values.split('-'):
                for nb in neigh_dict[int(added)]:
                    if str(nb) not in values.split('-'):
                        list_val = values.split('-')
                        sorted_list = [int(index) for index in list_val]
                        sorted_list.append(nb)
                        sorted_list.sort()
                        sorted_list_str = [str(val) for val in sorted_list]
                        new = '-'.join(sorted_list_str)
                        if sorted_list not in sorted_graphs:
                            subgraphs[len_subgraph].append(new)
                            sorted_graphs.append(sorted_list)
        if len(subgraphs[len_subgraph]) == 1:
            break
    return subgraphs