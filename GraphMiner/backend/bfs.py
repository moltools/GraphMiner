#!/usr/bin/env python3

def breadth_fs2(neigh_dict: dict, nodelist: list):
    '''
    Find all subgraphs that are present in a molecule

    input:
    nodelist - list containing the indeces (int) in the string of the atom
    neigh_dict - dictionary containing as key(int) the atom map numbers of the
    start atom and as value (list of int) the atom map numbers of the
    neighbouring atoms

    returns:
    subgraphs - dictionary containg as key (int) the subgraph length and as value
    (list of sets) with the subgraphs as sets
    '''
    len_subgraph = 1
    subgraphs = {}
    subgraphs[len_subgraph] = []
    # Generate the first list of keys
    for key in neigh_dict:
        subgraphs[len_subgraph].append({key})
    # Add in all unique regions for each length of subgraphs
    for index in range(len(nodelist)):
        len_subgraph += 1
        subgraphs[len_subgraph] = []
        sorted_graphs = []
        for valueset in subgraphs[len_subgraph - 1]:
            for atommapnum in valueset:
                for nb in neigh_dict[atommapnum]:
                    if nb not in valueset:
                        new = valueset.copy()
                        new.add(nb)
                        if new not in sorted_graphs:
                            subgraphs[len_subgraph].append(new)
                            sorted_graphs.append(new)
        if len(subgraphs[len_subgraph]) <= 1:
            break
    return subgraphs