#!/usr/bin/env python3

def list_of_nodes(insmiles:str, atom_list:list):
    '''
    Create a list of all nodes based on SMILES

    input:
    insmiles - SMILES format of a molecule (str)
    atom_list - list of all atoms (and merges) to check for in molecule]

    returns:
    node_list - list containing the indeces (int) in the string containing an atom
    '''
    node_list = []
    for index in range(len(insmiles)):
        if insmiles[index:index+2] in atom_list or insmiles[index] in atom_list:
            node_list.append(index)
    return node_list

def dict_of_nodes(insmiles:str, atom_list:list):
    '''
    Create a dictionary with neighboring SMILES
    
    input:
    insmiles - str of SMILES of a molecule
    atom_list - list of all atoms (and merges) to check for in molecule)
    
    returns:
    node_dict - dictionary containing as key (int) the index of the atom and as value
    a list of indeces (int) of the neighboring atoms.
    '''
    node_dict = {}
    for index in range(len(insmiles)):
        if insmiles[index:index+2] in atom_list and index+2 <= len(insmiles)-1:
            print('ADJUST FOR 2 LETTER ATOMS')
            if insmiles[index+2] == ')':
                continue
            elif index+2 <= len(insmiles)-1:
                node_dict[index] = [index+2]
        #Code for 1 letter atoms
        elif insmiles[index] in atom_list and index+1 <= len(insmiles)-1:
            #If the atom is followed by closing brackets, it does not bind to another
            #atom, making it unnecessary to add as a key value
            if insmiles[index+1] == ')':
                continue
            #If the atom is followed by opening brackets, it binds two atoms (at least),
            #making the group in between brackets and after brackets the two bound groups
            #How it is written now, does not work for lactose for example :(
            elif insmiles[index+1] == '(':
                if insmiles[index+2] == '=':
                    node_dict[index] = [index+3]
                    if insmiles[index+4] == '(' or ')' and index+5 <= len(insmiles)-1:
                        if insmiles[index+5] != '=' or '#':
                            node_dict[index].append(index+5)
                        else:
                            node_dict[index].append(index+5) 
                #For the elif part, same as above
                elif insmiles[index+2] != '=':
                    node_dict[index] = [index+2]
                    if insmiles[index+3] == '(' or ')' and index+4 <= len(insmiles)-1:
                        if insmiles[index+4] != '=' or '#':
                            node_dict[index].append(index+4)
                        else:
                            node_dict[index].append(index+4)
            #If there are no brackets, but just two atoms behind each other, 
            #they bind to each other, but to nothing else
            elif index+1 <= len(insmiles)-1:
                node_dict[index] = [index+1]
    return node_dict
        

def breadth_first_search(node_list:list, molsmiles:str, node_dict:dict):
    '''
    Perform Breadth First Search
    
    input:
    node_list - list containing the indeces (int) in the string containing an atom
    molsmiles - SMILES format of a molecule (str)
    node_dict - dictionary containing as key (int) the index of the atom and as value
    a list of indeces (int) of the neighboring atoms

    returns:
    vis - dictionary containing as key the index (int) of the atom and as value the 
    distance (int) to the first atom
    '''
    start_node = node_list[0]
    unvis = node_list
    curr_dist = 0
    queue = [start_node]
    distances = [0]
    vis = {}
    while len(queue) != 0:
        curr_node = queue[0]
        curr_dist = distances[0]
        vis.update({(curr_node):str(curr_dist)})
        if curr_node not in node_dict.keys():
            pass
        else:
            for node in node_dict[curr_node]:
                if node in queue:
                    continue
                elif (node) not in unvis:
                    continue
                else:
                    queue.append(node)
                    distances.append(int(curr_dist)+1)
                    unvis.remove(((node)))
        queue.pop(0)
        distances.pop(0)
    node_list.pop(0)
    for nodes in unvis:
        vis.update({nodes:'-1'})
    return vis

    
