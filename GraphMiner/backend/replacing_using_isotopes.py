#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem
from rdkit.Chem import AllChem
import string

#It does not retrn the correct ones :(
#Don't know how/what/when

def combining(molsmiles:str):
    moll = Chem.MolFromSmiles(molsmiles)
    if moll.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')):
        indeces = moll.GetSubstructMatches(Chem.MolFromSmiles('C(=O)O'))
        for number in range(len(indeces)):
            patt = Chem.MolFromSmiles('C(=O)O')
            repl = Chem.MolFromSmiles('C[1*]')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    if moll.HasSubstructMatch(Chem.MolFromSmiles('CO')):
        indeces = moll.GetSubstructMatches(Chem.MolFromSmiles('CO'))
        for number in range(len(indeces)):
            patt = Chem.MolFromSmiles('CO')
            repl = Chem.MolFromSmiles('C[2*]')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    mts = Chem.MolToSmiles(moll)
    return mts, indeces

def returning(mts:str, indeces):
    moll = Chem.MolFromSmiles(mts)
    if moll.HasSubstructMatch(Chem.MolFromSmiles('C[1*]')):
        for number in range(len(indeces)):
            patt = Chem.MolFromSmiles('C[1*]')
            repl = Chem.MolFromSmiles('C(=O)O')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    if moll.HasSubstructMatch(Chem.MolFromSmiles('C[2*]')):
        for number in range(len(indeces)):
            patt = Chem.MolFromSmiles('C[2*]')
            repl = Chem.MolFromSmiles('CO')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    new_smile = Chem.MolToSmiles(moll)
    return new_smile