#!/usr/bin/env python3

#IMPORTS
import rdkit
from rdkit import Chem


# def prep_replacing(insmiles):
#     inmol = Chem.MolFromSmiles(insmiles)
#     return inmol


# def replacing_COO(moll, replacements:dict):
#     list_of_C = [Cidx[0] for Cidx in
#                  moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
#     rwmol = rdkit.Chem.rdchem.RWMol()
#     rwmol.InsertMol(moll)
#     COOlist = []
#     query = Chem.MolFromSmiles('C(=O)O')
#     idx_COO = moll.GetSubstructMatches(query)
#     for idx__COO in idx_COO:
#         sort_idx = sorted(idx__COO)
#         for indiv_idx in idx__COO:
#             if indiv_idx not in list_of_C:
#                 COOlist.append(indiv_idx)
#         replacements[sort_idx[0]] = sort_idx[1:]
#     # Actual Replacement
#     COOlist.sort()
#     rem_atoms = 0
#     for number in COOlist:
#         adj_num = number-rem_atoms
#         rwmol.RemoveAtom(adj_num)
#         if '.' in Chem.MolToSmiles(rwmol):
#             rwmol.RemoveBond(adj_num-1, adj_num)
#             rwmol.AddBond(adj_num - 1, adj_num, rdkit.Chem.rdchem.BondType.SINGLE)
#         rem_atoms += 1
#     inputsmls = Chem.MolToSmiles(rwmol)
#     moll = Chem.MolFromSmiles(inputsmls)
#     return moll, replacements


# def replacing_CO(moll, replacements:dict):
#     list_of_C = [Cidx[0] for Cidx in
#                  moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
#     rwmol = rdkit.Chem.rdchem.RWMol()
#     rwmol.InsertMol(moll)
#     COlist = []
#     idx_CO = moll.GetSubstructMatches(Chem.MolFromSmiles('CO'))
#     for idx__CO in idx_CO:
#         sort_idx = sorted(idx__CO)
#         for indiv_idx in idx__CO:
#             if indiv_idx not in list_of_C:
#                 COlist.append(indiv_idx)
#         replacements[sort_idx[0]] = sort_idx[1:]
#     # Actual Replacement
#     COlist.sort()
#     rem_atoms = 0
#     for number in COlist:
#         adj_num = number-rem_atoms
#         rwmol.RemoveAtom(adj_num)
#         if '.' in Chem.MolToSmiles(rwmol):
#             rwmol.RemoveBond(adj_num-1, adj_num)
#             rwmol.AddBond(adj_num - 1, adj_num, rdkit.Chem.rdchem.BondType.SINGLE)
#         rem_atoms += 1
#     inputsmls = Chem.MolToSmiles(rwmol)
#     moll = Chem.MolFromSmiles(inputsmls)
#     return moll, replacements


# def replacing_C_O(moll, replacements:dict):
#     list_of_C = [Cidx[0] for Cidx in
#                  moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
#     rwmol = rdkit.Chem.rdchem.RWMol()
#     rwmol.InsertMol(moll)
#     COlist = []
#     idx_CO = moll.GetSubstructMatches(Chem.MolFromSmiles('C=O'))
#     for idx__CO in idx_CO:
#         sort_idx = sorted(idx__CO)
#         for indiv_idx in idx__CO:
#             if indiv_idx not in list_of_C:
#                 COlist.append(indiv_idx)
#         replacements[sort_idx[0]] = sort_idx[1:]
#     # Actual Replacement
#     COlist.sort()
#     rem_atoms = 0
#     for number in COlist:
#         adj_num = number-rem_atoms
#         rwmol.RemoveAtom(adj_num)
#         if '.' in Chem.MolToSmiles(rwmol):
#             rwmol.RemoveBond(adj_num-1, adj_num)
#             rwmol.AddBond(adj_num - 1, adj_num, rdkit.Chem.rdchem.BondType.SINGLE)
#         rem_atoms += 1
#     inputsmls = Chem.MolToSmiles(rwmol)
#     moll = Chem.MolFromSmiles(inputsmls)
#     return moll, replacements
