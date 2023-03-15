#!/usr/bin/env python3

#IMPORTS
import rdkit
from rdkit import Chem


def starting_rwmol(inputsmls):
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(Chem.MolFromSmiles(inputsmls))
    moll = Chem.MolFromSmiles(inputsmls)
    replacements = {}
    list_of_C = [Cidx[0] for Cidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
    # if rwmol.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
    #     idx_COO = rwmol.GetSubstructMatches(Chem.MolFromSmiles('C(=O)O'))
    #     for idx__COO in idx_COO:
    #         new_list = []
    #         for indiv_idx in idx__COO:
    #             if indiv_idx in list_of_C:
    #                 start_idx = indiv_idx
    #             elif indiv_idx not in list_of_C:
    #                 new_list.append(indiv_idx)
    #         replacements[start_idx] = new_list
    #     # Actual Replacement
    #     for number in range(len(idx_COO)):
    #         patt = Chem.MolFromSmiles('C(=O)O')
    #         repl = Chem.MolFromSmiles('C')
    #         repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
    #         moll = repl_str[0]
    # mts = Chem.MolToSmiles(moll)
    return mts, replacements