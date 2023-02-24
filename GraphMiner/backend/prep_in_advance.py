#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem
from rdkit.Chem import AllChem

def select_on_size(smile_mol):
    ''' 
    Only return SMILES if number of heavy atoms below 20
    
    input: 
    smile_mol - SMILES format of a molecule (str)
    
    returns:
    smile_mol - SMILES format of a molecule (str)'''
    mol1 = Chem.MolFromSmiles(smile_mol)
    heavy_atoms = mol1.GetNumHeavyAtoms()
    if heavy_atoms <= 20: #Change to 40
        return smile_mol


def combine_basic_substructures(molsmiles):
    print(molsmiles)
    moll = Chem.MolFromSmiles(molsmiles)
    if 'C(=O)O' in molsmiles:
        patt = Chem.MolFromSmiles('C(=O)O')
        repl = Chem.MolFromSmiles('C')
        for atom in atomIter:
            atom.GetIdx()
        repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
        mts = Chem.MolToSmiles(repl_str[0])
    else:
        mts = molsmiles
    return mts