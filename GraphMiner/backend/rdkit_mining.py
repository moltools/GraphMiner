#!/usr/bin/env python3

from rdkit import Chem

def subgraph_miner(Smiles:str):
    '''
    Retrieve the subgraphs from a molecule in smiles format
    
    input:
    Smiles - molecule in canonical SMILES format
    
    returns:
    subgraphs - tuple of lists containing tuples. Tuples indicate
    the number of atoms in the subgraph, each list has tuples of the same length.
    m1 - molecule in mol format from RDKit'''
    m1 = Chem.MolFromSmiles(Smiles)
    heavy_atoms = m1.GetNumHeavyAtoms()
    subgraphs = Chem.FindAllSubgraphsOfLengthMToN(m1, 1, m1.GetNumHeavyAtoms()-1)
    return subgraphs, m1, heavy_atoms

def sub_to_smiles(subgraphs:tuple, mol1):
    '''
    Turn subgraphs retrieved into SMILES
    
    input:
    subgraphs - tuple of lists containing tuples. Tuples indicate
    the number of atoms in the subgraph, each list has tuples of the same length.
    mol1 - molecule in mol format from RDKit

    returns:
    list_of_smiles - list containing all subgraphs in SMILES format
    '''
    list_of_smiles = []
    for size in subgraphs:
        for atoms in size:
            sub_smiles = Chem.MolFragmentToSmiles(mol1, atoms)
            print(sub_smiles)
            list_of_smiles.append(sub_smiles)
    return list_of_smiles