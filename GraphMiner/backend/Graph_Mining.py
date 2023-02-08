#!/usr/bin/env python3

def list_of_nodes(insmiles:str, atom_list:list):
    '''
    Create a list of all nodes based on SMILES

    input:
    insmiles - SMILES format of a molecule
    atom_list - list of all atoms (and merges) to check for in molecule
    '''
    node_list = []
    for index in range(len(insmiles)):
        if insmiles[index:index+2] in atom_list:
            node_list.append(insmiles[index:index+2])
        elif insmiles[index] in atom_list:
            node_list.append(insmiles[index])
    return node_list

def dict_of_nodes(insmiles:str, atom_list:list, node_list:list):
    node_dict = {}
    # for index in range(len(insmiles)):
    #     if insmiles[index:index+2] in atom_list:
    #         node_dict[]
    #     elif insmiles[index] in atom_list:
    #         node_dict[]


    return

        

def breadth_first_search(node_list:list, molsmiles:str):
    start_node = node_list[0]
    unvis = node_list
    curr_dist = 0
    queue = [start_node]
    distances = [0]
    vis = {}
    # while len(queue) != 0:
    #     curr_node = queue[0]
    #     curr_dist = distances[0]
    #     vis.update({(curr_node):str(curr_dist)})
    #     if curr_node not in dict_node.keys():
    #         pass
    #     else:
    #         for node in dict_node[curr_node]:
    #             if node in queue:
    #                 continue
    #             elif (node) not in unvis:
    #                 continue
    #             else:
    #                 queue.append(node)
    #                 distances.append(int(curr_dist)+1)
    #                 unvis.remove(((node)))
    #     queue.pop(0)
    #     distances.pop(0)
    # node_list.pop(0)
    # for nodes in unvis:
    #     vis.update({nodes:'-1'})
    return vis

    
