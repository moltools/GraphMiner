#!/usr/bin/env python3

#IMPORTS
import rdkit
from rdkit import Chem

def set_atommapnum(molecule):
    '''
    Give each atom a unique AtomMapNumber

    input:
    molecule - molecule in MOL format from RDKit

    output:
    None, but atommap numbers are stored to the inidivual atoms in the molecule
    '''
    atommapnum = 0
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(atommapnum)
        atommapnum += 1
    return

def repl_atommap_COO(moll, replacements:dict):
    '''
    Replace the substructure C(=O)O by C using RWMol and AtomMapNumbers

    input:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule

    returns:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    -- but both are adjusted with the replaced C(=O)O substructures
    '''
    list_of_C = [Cidx[0] for Cidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    COOlist = []
    query = Chem.MolFromSmiles('C(=O)O')
    idx_COO = moll.GetSubstructMatches(query)
    for idx__COO in idx_COO:
        new_list = []
        for indiv_idx in idx__COO:
            if indiv_idx not in list_of_C:
                if indiv_idx not in COOlist:
                    COOlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_C:
                C_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[C_atom] = new_list
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

def repl_atommap_NCO(moll, replacements:dict):
    '''
    Replace the substructure NC=O by C using RWMol and AtomMapNumbers

    input:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule

    returns:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    -- but both are adjusted with the replaced NC=O substructures
    '''
    list_of_C = [Cidx[0] for Cidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    NCOlist = []
    query = Chem.MolFromSmiles('NC(=O)')
    idx_NCO = moll.GetSubstructMatches(query)
    for idx__NCO in idx_NCO:
        new_list = []
        for indiv_idx in idx__NCO:
            if indiv_idx not in list_of_C:
                if indiv_idx not in NCOlist:
                    NCOlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_C:
                C_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[C_atom] = new_list
    # Actual Replacement
    NCOlist.sort()
    rem_atoms = 0
    for number in NCOlist:
        adj_num = number-rem_atoms
        rwmol.RemoveAtom(adj_num)
        if '.' in Chem.MolToSmiles(rwmol):
            rwmol.RemoveBond(adj_num-1, adj_num)
            rwmol.AddBond(adj_num - 1, adj_num, rdkit.Chem.rdchem.BondType.SINGLE)
        rem_atoms += 1
    inputsmls = Chem.MolToSmiles(rwmol)
    moll = Chem.MolFromSmiles(inputsmls)
    return moll, replacements

def repl_atommap_CO(moll, replacements:dict):
    '''
    Replace the substructure CO by C using RWMol and AtomMapNumbers

    input:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule

    returns:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    -- but both are adjusted with the replaced CO substructures
    '''
    list_of_C = [Cidx[0] for Cidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    COlist = []
    query = Chem.MolFromSmiles('CO')
    idx_CO = moll.GetSubstructMatches(query)
    for idx__CO in idx_CO:
        new_list = []
        for indiv_idx in idx__CO:
            if indiv_idx not in list_of_C:
                if indiv_idx not in COlist:
                    COlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_C:
                C_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[C_atom] = new_list
    # Actual Replacement
    COlist.sort()
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

def repl_atommap_C_O(moll, replacements:dict):
    '''
    Replace the substructure C=O by C using RWMol and AtomMapNumbers

    input:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule

    returns:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    -- but both are adjusted with the replaced C=O substructures
    '''
    list_of_C = [Cidx[0] for Cidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('C'))]
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    COlist = []
    query = Chem.MolFromSmiles('C=O')
    idx_CO = moll.GetSubstructMatches(query)
    for idx__CO in idx_CO:
        new_list = []
        for indiv_idx in idx__CO:
            if indiv_idx not in list_of_C:
                if indiv_idx not in COlist:
                    COlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_C:
                C_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[C_atom] = new_list
    # Actual Replacement
    COlist.sort()
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

def return_replaced(repl_dicts:dict, index_dicts:dict):
    '''
    Return the AtomMapNumbers of the replaced substructure in found substructures

    input:
    repl_dicts - dictionary with as key the AtomMapNumber (int) of the C that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    index_dicts - dictionary with as key the length (int) of the subgraph and
    as value a list of all subgraphs (str) as AtomMapNumbers divided by dashes

    returns:
    index_dicts - dictionary with as key the length (int) of the subgraph and
    as value a list of all subgraphs (str) as AtomMapNumbers divided by dashes
    with the replaced substructure AtomMapNumbers
    '''
    for subgraphlength in index_dicts:
        for value in repl_dicts:
            for subgraph in index_dicts[subgraphlength]:
                if str(value) not in subgraph:
                    continue
                index = index_dicts[subgraphlength].index(subgraph)
                sorted_list = [int(idx) for idx in subgraph.split('-')]
                if str(value) in subgraph.split('-'):
                    sorted_list += (repl_dicts[value])
                    sorted_list.sort()
                sorted_list_str = [str(val) for val in sorted_list]
                r_new = '-'.join(sorted_list_str)
                index_dicts[subgraphlength][index] = r_new
    return index_dicts
