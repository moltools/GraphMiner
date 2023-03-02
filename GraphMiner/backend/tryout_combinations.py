#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem
from rdkit.Chem import AllChem
import string

#It does not retrn the correct ones :(
#Don't know how/what/when

def combining(molsmiles:str):
    print(' ')
    print(molsmiles)
    moll = Chem.MolFromSmiles(molsmiles)
    indeces = moll.GetSubstructMatches(Chem.MolFromSmiles('CO'))
    for number in range(len(indeces)):
        patt = Chem.MolFromSmiles('CO')
        repl = Chem.MolFromSmiles('Br')
        repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
        moll = repl_str[0]
    mts = Chem.MolToSmiles(moll)
    return mts, indeces, moll

def returning(mts:str, indeces, moll):
    for number in range(len(indeces)):
        patt = Chem.MolFromSmiles('Br')
        repl = Chem.MolFromSmiles('CO')
        repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
        moll = repl_str[0]
    new_smile = Chem.MolToSmiles(moll)
    return new_smile