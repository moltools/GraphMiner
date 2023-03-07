#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem
from rdkit.Chem import AllChem
import string

#It does not retrn the correct ones :(
#Don't know how/what/when

def combining(molsmiles:str):
    '''
    Replace substructures by an isotope of C

    input:
    molsmiles - SMILES format of a molecule (str)

    returns:
    mts - SMILES format of a molecule, with combined substructures (str)
    indeces - list of tuples containing tuples of the indices of the substructures
    '''
    moll = Chem.MolFromSmiles(molsmiles)
    indeces = []
    if moll.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')):
        indeces.append(moll.GetSubstructMatches(Chem.MolFromSmiles('C(=O)O')))
        for number in range(len(indeces)):
            patt = Chem.MolFromSmiles('C(=O)O')
            repl = Chem.MolFromSmiles('C[1*]')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    if moll.HasSubstructMatch(Chem.MolFromSmiles('CO')):
        indeces.append(moll.GetSubstructMatches(Chem.MolFromSmiles('CO')))
        for number in range(len(indeces)):
            patt = Chem.MolFromSmiles('CO')
            repl = Chem.MolFromSmiles('C[2*]')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    mts = Chem.MolToSmiles(moll)
    return mts, indeces

def returning(mts:str, indeces:list):
    '''
    Replace isotope of C by the replaced substructure

    input:
    mts - SMILES format of a molecule, with combined substructures (str)
    indeces - list of tuples containing tuples of the indices of the substructures

    returns:
    new_smile - SMILES format of a molecule with backreplaced substructures
    '''
    moll = Chem.MolFromSmiles(mts)
    if moll.HasSubstructMatch(Chem.MolFromSmiles('C[1*]')):
        for number in range(len(indece[1])):
            patt = Chem.MolFromSmiles('C[1*]')
            repl = Chem.MolFromSmiles('C(=O)O')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    if moll.HasSubstructMatch(Chem.MolFromSmiles('C[2*]')):
        for number in range(len(indeces[2])):
            patt = Chem.MolFromSmiles('C[2*]')
            repl = Chem.MolFromSmiles('CO')
            repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
            moll = repl_str[0]
    new_smile = Chem.MolToSmiles(moll)
    return new_smile