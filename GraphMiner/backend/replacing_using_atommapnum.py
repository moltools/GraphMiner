#!/usr/bin/env python3

#IMPORTS
import rdkit
from rdkit import Chem

def set_atommapnum(molecule):
    atommapnum = 0
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(atommapnum)
        atommapnum += 1
    return

def repl_atommap_COO(moll, replacements:dict):
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
    print(replacements)
    return moll, replacements

def return_replaced(repl_dicts:dict, index_dicts:dict):
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
