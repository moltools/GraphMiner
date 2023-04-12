#!/usr/bin/env python3

def breadth_fs(nodelist:list, neigh_dict:dict):
    '''
    Find all subgraphs that are present
    
    input:
    nodelist - list containing the indeces (int) in the string of the atom
    neigh_dict - dictionary containing as key(int) the index of the start atom
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
    # print('step1')
    #Add in all unqiue regions for each length of subgraphs
    for index in range(len(nodelist)):
        len_subgraph += 1
        subgraphs[len_subgraph] = []
        sorted_graphs = []
        # print(index)
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
        if len(subgraphs[len_subgraph]) <= 1:
            break
    return subgraphs


def breadth_fs2(neigh_dict: dict, nodelist: list):
    '''
    Find all subgraphs that are present

    input:
    nodelist - list containing the indeces (int) in the string of the atom
    neigh_dict - dictionary containing as key(int) the index of the start atom
    and as value (list of int) the indices of the neighbouring atoms

    returns:
    subgraphs - dictionary containg as key (int) the subgraph length and as value
    (list of str) with the subgraphs as strings of indeces divided by '-'
    '''
    len_subgraph = 1
    subgraphs = {}
    subgraphs[len_subgraph] = []
    # Generate the first list of keys
    for key in neigh_dict:
        subgraphs[len_subgraph].append({key})
    # Add in all unqiue regions for each length of subgraphs
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

# def breadth_fs2(molecule, dict_node, list_node):
#     '''xx
#     '''
#     print(list_node)
#     for atom in molecule.GetAtoms():
#         start_node = atom.GetAtomMapNum()
#         unvis = list_node
#         curr_dist = 0
#         queue = [start_node]
#         distances = [0]
#         vis = {}
#         while len(queue) != 0:
#             print(queue)
#             curr_node = queue[0]
#             curr_dist = distances[0]
#             vis.update({int(curr_node):str(curr_dist)})
#             if curr_node not in dict_node.keys():
#                 pass
#             else:
#                 for node in dict_node[curr_node]:
#                     if node in queue:
#                         continue
#                     # elif int(node) not in unvis:
#                     #     continue
#                     else:
#                         queue.append(node)
#                         distances.append(int(curr_dist)+1)
#                         # unvis.remove((int(node)))
#             queue.pop(0)
#             distances.pop(0)
#         # print(vis)
#     return vis