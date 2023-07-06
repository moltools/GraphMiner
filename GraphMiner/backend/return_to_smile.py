#!/usr/bin/env python3

from rdkit import Chem

def rdkit_smiles2(sub_graphs:dict, smilesmol, zeromol):
    '''
    Convert the indices of subgraphs to SMILES subgraphs

    input:
    subgraphs - dictionary containing as key (int) the subgraph length and as value
    (set of int) with the subgraphs as integers
    smilesmol - MOL format of a molecule (mol)

    returns:
    mol_graphs - dictionary containing as key (int) the subgraph length and as value
    (list of str) the subgraphs as strings in SMILES format
    '''
    mol_graphs = {}
    smiles_graphs = {}
    atommapdict = {}
    dot_structure = 0
    for atom in smilesmol.GetAtoms():
        atommapdict[atom.GetAtomMapNum()] = atom.GetIdx()
        atom.SetAtomMapNum(0)
    for subgraph_length in sub_graphs:
        smiles_graphs[subgraph_length] = []
        for subgraphset in sub_graphs[subgraph_length]:
            for atommapnum in subgraphset:
                if atommapnum in atommapdict.keys():
                    subgraphset.remove(atommapnum)
                    subgraphset.add(atommapdict[atommapnum])
            subgraphlist = list(subgraphset)
            if '.' in Chem.MolFragmentToSmiles(smilesmol, subgraphlist):
                dot_structure += 1
                continue
            smiles_graphs[subgraph_length].append(
                Chem.MolFragmentToSmiles(zeromol, subgraphlist))
    print('dot_structure', dot_structure)
    return smiles_graphs, mol_graphs

def rdkit_smiles3(sub_graphs:dict, smilesmol):
    '''
    Convert the indices of subgraphs to SMILES subgraphs

    input:
    subgraphs - dictionary containing as key (int) the subgraph length and as value
    (set of int) with the subgraphs as integers
    smilesmol - MOL format of a molecule (mol)

    returns:
    mol_graphs - dictionary containing as key (int) the subgraph length and as value
    (list of str) the subgraphs as strings in SMILES format
    '''
    mol_graphs = []
    atommapdict = {}
    for atom in smilesmol.GetAtoms():
        if atom.GetAtomMapNum() != atom.GetIdx():
            atommapdict[atom.GetAtomMapNum()] = atom.GetIdx()
    if len(atommapdict) >0:
        for subgraphset in sub_graphs:
            for atommapnum in subgraphset:
                if atommapnum in atommapdict.keys():
                    subgraphset.remove(atommapnum)
                    subgraphset.add(atommapdict[atommapnum])
                subgraphlist = list(subgraphset)
            mol_graphs.append(Chem.MolFragmentToSmiles(smilesmol, subgraphlist))
    else:
        for subgraphset in sub_graphs:
            subgraphlist = list(subgraphset)
            mol_graphs.append(Chem.MolFragmentToSmiles(smilesmol, subgraphlist))
    return mol_graphs

