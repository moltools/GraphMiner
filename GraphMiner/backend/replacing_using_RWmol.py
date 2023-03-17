#!/usr/bin/env python3

#IMPORTS
import rdkit
from rdkit import Chem


def prep_replacing(insmiles):
    inmol = Chem.MolFromSmiles(insmiles)
    Clist = [Cidx[0] for Cidx in
                 inmol.GetSubstructMatches(Chem.MolFromSmiles('C'))]
    return Clist, inmol


def replacing_COO(moll, list_of_C:list, replacements:dict):
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    COOlist = []
    idx_COO = moll.GetSubstructMatches(Chem.MolFromSmiles('C(=O)O'))
    for idx__COO in idx_COO:
        new_list = []
        for indiv_idx in idx__COO:
            if indiv_idx in list_of_C:
                start_idx = indiv_idx
            elif indiv_idx not in list_of_C:
                new_list.append(indiv_idx)
                COOlist.append(indiv_idx)
        replacements[start_idx] = new_list
    # Actual Replacement
    COOlist.sort()
    rem_atoms = 0
    for number in COOlist:
        adj_num = number-rem_atoms
        rwmol.RemoveAtom(adj_num)
        if '.' in Chem.MolToSmiles(rwmol):
            rwmol.RemoveBond(adj_num-1, adj_num)
            rwmol.AddBond(adj_num - 1, adj_num, rdkit.Chem.rdchem.BondType.SINGLE)
        rem_atoms += 1
    inputsmls = Chem.MolToSmiles(rwmol)
    moll = Chem.MolFromSmiles(inputsmls)
    return moll, replacements


def replacing_CO(moll, list_of_C:list, replacements:dict):
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    COlist = []
    idx_CO = moll.GetSubstructMatches(Chem.MolFromSmiles('CO'))
    for idx__CO in idx_CO:
        new_list = []
        for indiv_idx in idx__CO:
            if indiv_idx in list_of_C:
                start_idx = indiv_idx
            elif indiv_idx not in list_of_C:
                new_list.append(indiv_idx)
                COlist.append(indiv_idx)
        replacements[start_idx] = new_list
    # Actual Replacement
    COlist.sort()
    print(COlist)
    rem_atoms = 0
    for number in COlist:
        adj_num = number-rem_atoms
        rwmol.RemoveAtom(adj_num)
        print(Chem.MolToSmiles(rwmol))
        if '.' in Chem.MolToSmiles(rwmol):
            print('dot found')
            print(rwmol.GetBondBetweenAtoms(adj_num-1, adj_num))
            rwmol.RemoveBond(adj_num-1, adj_num)
            rwmol.AddBond(adj_num - 1, adj_num, rdkit.Chem.rdchem.BondType.SINGLE)
        rem_atoms += 1
    inputsmls = Chem.MolToSmiles(rwmol)
    moll = Chem.MolFromSmiles(inputsmls)
    return moll, replacements


def replacing_C_O(moll, list_of_C:list, replacements:dict):
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    COlist = []
    idx_CO = moll.GetSubstructMatches(Chem.MolFromSmiles('C=O'))
    for idx__CO in idx_CO:
        new_list = []
        for indiv_idx in idx__CO:
            if indiv_idx in list_of_C:
                start_idx = indiv_idx
            elif indiv_idx not in list_of_C:
                new_list.append(indiv_idx)
                COlist.append(indiv_idx)
        replacements[start_idx] = new_list
    # Actual Replacement
    COlist.sort()
    print(COlist)
    rem_atoms = 0
    for number in COlist:
        adj_num = number-rem_atoms
        rwmol.RemoveAtom(adj_num)
        if '.' in Chem.MolToSmiles(rwmol):
            rwmol.RemoveBond(adj_num-1, adj_num)
            rwmol.AddBond(adj_num - 1, adj_num, rdkit.Chem.rdchem.BondType.SINGLE)
        rem_atoms += 1
    inputsmls = Chem.MolToSmiles(rwmol)
    moll = Chem.MolFromSmiles(inputsmls)
    return moll, replacements
