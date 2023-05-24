#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem

def select_on_size(smile_mol:str, sizelimit:int):
    ''' 
    Only return SMILES if number of heavy atoms below 20
    
    input: 
    smile_mol - SMILES format of a molecule (str)
    
    returns:
    smile_mol - SMILES format of a molecule (str)
    '''
    mol1 = Chem.MolFromSmiles(smile_mol)
    heavy_atoms = mol1.GetNumHeavyAtoms()
    if heavy_atoms <= sizelimit: #Change to 40
        return smile_mol

def select_mol(molsmile: str):
    if molsmile == None:
        return
    if '.' in molsmile:
        selected = molsmile.split('.')
        if len(selected) == 0:
            return
        elif len(selected) == 1:
            return Chem.MolFromSmiles(selected[0])
        elif len(selected) > 1:
            return selected
    return Chem.MolFromSmiles(molsmile)
