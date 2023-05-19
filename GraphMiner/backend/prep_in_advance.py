#!/usr/bin/env python3

#IMPORT STATEMENTS
from rdkit import Chem
from rdkit.Chem import AllChem
import string

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
        splitmol = molsmile.split('.')
        selected = []
        for part in splitmol:
            if Chem.MolFromSmiles(part).GetNumHeavyAtoms() > 5:
                selected.append(part)
        if len(selected) == 0:
            return
        elif len(selected) == 1:
            return Chem.MolFromSmiles(selected[0])
        elif len(selected) > 1:
            return selected
    return Chem.MolFromSmiles(molsmile)

def remove_atom_charges(mol: Chem.Mol) -> Chem.Mol:
    """
    Removes charges from molecule.
    Arguments
    ---------
    target : Chem.Mol
        RDKit `Chem.Mol` object to query for patterns.
    Returns
    -------
    target : Chem.Mol
        RDKit `Chem.Mol` object to query for patterns with neutralized charges.
    Source: rdkit.org/docs/Cookbook.html
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            if atom.GetIsAromatic() or any([a.GetIsAromatic() for a in atom.GetNeighbors()]):
                # Skip removing charges for aromatic-adjacent atoms which causes valence issues.
                continue
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def ed_remove_atom_charges(mol: Chem.Mol) -> Chem.Mol:
    """
    Removes charges from molecule.
    Arguments
    ---------
    target : Chem.Mol
        RDKit `Chem.Mol` object to query for patterns.
    Returns
    -------
    target : Chem.Mol
        RDKit `Chem.Mol` object to query for patterns with neutralized charges.
    Source: rdkit.org/docs/Cookbook.html
    """
    pattern = Chem.MolFromSmarts("[+]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        print('USED')
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            if atom.GetIsAromatic() or any([a.GetIsAromatic() for a in atom.GetNeighbors()]):
                # Skip removing charges for aromatic-adjacent atoms which causes valence issues.
                continue
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            print(hcount)
            atom.SetFormalCharge(0)
            if hcount - chg > 0:
                atom.SetNumExplicitHs(hcount - chg)
            else:
                atom.SetNumExplicitHs(0)
            atom.UpdatePropertyCache()
    return mol


# def combine_basic_substructures(molsmiles:str):
#     '''
#     Replace substructures by C and store replacements
#
#     input:
#     molsmiles - SMILES format of a molecule (str)
#
#     returns:
#     mts - SMILES format of a molecule, with combined substructures (str)
#     replacements - dictionary containing as key (int) the index of the remaining
#     atom of the substructure and as value (list of int) the indeces of the replaced atoms
#     '''
#     moll = Chem.MolFromSmiles(molsmiles)
#     replacements = {}
#     list_of_C = [Cidx[0] for Cidx in moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
#     if moll.HasSubstructMatch(Chem.MolFromSmiles('C(=O)O')) == True:
#         idx_COO = moll.GetSubstructMatches(Chem.MolFromSmiles('C(=O)O'))
#         for idx__COO in idx_COO:
#             new_list = []
#             for indiv_idx in idx__COO:
#                 if indiv_idx in list_of_C:
#                     start_idx = indiv_idx
#                 elif indiv_idx not in list_of_C:
#                     new_list.append(indiv_idx)
#             replacements[start_idx] = new_list
#         #Actual Replacement
#         for number in range(len(idx_COO)):
#             patt = Chem.MolFromSmiles('C(=O)O')
#             repl = Chem.MolFromSmiles('C')
#             repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
#             moll = repl_str[0]
#     if moll.HasSubstructMatch(Chem.MolFromSmiles('CO')) == True:
#         idx_CO = moll.GetSubstructMatches(Chem.MolFromSmiles('CO'))
#         for idx__CO in idx_CO:
#             new_list = []
#             for indiv_idx in idx__CO:
#                 if indiv_idx in list_of_C:
#                     start_idx = indiv_idx
#                 elif indiv_idx not in list_of_C:
#                     new_list.append(indiv_idx)
#             if start_idx not in replacements.keys():
#                 replacements[start_idx] = new_list
#         for number in range(len(idx_CO)):
#             patt = Chem.MolFromSmiles('CO')
#             repl = Chem.MolFromSmiles('C')
#             repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
#             moll = repl_str[0]
#     if moll.HasSubstructMatch(Chem.MolFromSmiles('C(=O)')) == True:
#         idx_C_O = moll.GetSubstructMatches(Chem.MolFromSmiles('C(=O)'))
#         for idx__C_O in idx_C_O:
#             new_list = []
#             for indiv_idx in idx__C_O:
#                 if indiv_idx in list_of_C:
#                     start_idx = indiv_idx
#                 elif indiv_idx not in list_of_C:
#                     new_list.append(indiv_idx)
#             if start_idx not in replacements.keys():
#                 replacements[start_idx] = new_list
#         for number in range(len(idx_C_O)):
#             patt = Chem.MolFromSmiles('C(=O)')
#             repl = Chem.MolFromSmiles('C')
#             repl_str = AllChem.ReplaceSubstructs(moll, patt, repl)
#             moll = repl_str[0]
#     mts = Chem.MolToSmiles(moll)
#     return mts, replacements
#
#
# def return_basic_substructures(repl_dicts:dict, index_dicts:dict):
#     '''
#     Add in substrucures again
#
#     input:
#     repl_dicts - dictionary containing as key (int) the index of the remaining
#     atom of the substructure and as value (list of int) the indeces of the replaced atoms
#     index_dicts - dictionary containg as key (int) the subgraph length and as value
#     (list of str) with the subgraphs as strings of indeces divided by '-'
#
#     returns:
#     index_dicts - dictionary containg as key (int) the subgraph length and as value
#     (list of str) with the subgraphs as strings of indeces divided by '-', but with
#     substructures indeces placed back in again
#     '''
#     for subgraphlength in index_dicts:
#         for value in repl_dicts:
#             for subgraph in index_dicts[subgraphlength]:
#                 if subgraph not in index_dicts[subgraphlength]:
#                     continue
#                 index = index_dicts[subgraphlength].index(subgraph)
#                 sorted_list = [int(idx) for idx in subgraph.split('-')]
#                 # Replace all values that are higher than the one for which the
#                 # subgraph will be added in
#                 length_input = len(repl_dicts[value])
#                 for ind in range(len(sorted_list)):
#                     if sorted_list[ind] > value:
#                         sorted_list[ind] = sorted_list[ind] + length_input
#                     else:
#                         for val in repl_dicts[value]:
#                             if value > sorted_list[ind] >= val:
#                                 sorted_list[ind] = value
#
#                 sorted_list_str = [str(val) for val in sorted_list]
#                 new = '-'.join(sorted_list_str)
#                 index_dicts[subgraphlength][index] = new
#                 # print(index_dicts)
#                 #Add in the earlier combined subgraphs
#                 if str(value) in subgraph.split('-'):
#                     sorted_list += (repl_dicts[value])
#                     sorted_list.sort()
#                 sorted_list_str = [str(val) for val in sorted_list]
#                 r_new = '-'.join(sorted_list_str)
#                 index_dicts[subgraphlength][index] = r_new
#     return index_dicts