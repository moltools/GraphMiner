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
    print(idx_COO)
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
        adj_num = number - rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        # print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
    return moll, replacements

def repl_atommap_POOO(moll, replacements:dict):
    '''
    Replace the substructure P(=O)(O)O by P using RWMol and AtomMapNumbers

    input:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the P that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule

    returns:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the P that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    -- but both are adjusted with the replaced P(=O)(O)O substructures
    '''
    list_of_P = [Pidx[0] for Pidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('P'))]
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    POOOlist = []
    query = Chem.MolFromSmiles('P(=O)(O)O')
    idx_POOO = moll.GetSubstructMatches(query)
    for idx__POOO in idx_POOO:
        new_list = []
        for indiv_idx in idx__POOO:
            if indiv_idx not in list_of_P:
                if indiv_idx not in POOOlist:
                    POOOlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_P:
                P_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[P_atom] = new_list
    # Actual Replacement
    POOOlist.sort()
    rem_atoms = 0
    for number in POOOlist:
        adj_num = number - rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        # print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
    return moll, replacements

def repl_atommap_NOCO(moll, replacements:dict):
    '''
    Replace the substructure NO by C using RWMol and AtomMapNumbers

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
    NOCOlist = []
    query = Chem.MolFromSmiles('N(O)C(=O)')
    idx_NOCO = moll.GetSubstructMatches(query)
    for idx__NOCO in idx_NOCO:
        new_list = []
        for indiv_idx in idx__NOCO:
            if indiv_idx not in list_of_C:
                if indiv_idx not in NOCOlist:
                    NOCOlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_C:
                C_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[C_atom] = new_list
    # Actual Replacement
    NOCOlist.sort()
    rem_atoms = 0
    for number in NOCOlist:
        adj_num = number - rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        # print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
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
        adj_num = number - rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        # print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
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
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
                if indiv_idx not in COlist:
                    COlist.append(indiv_idx)
            elif indiv_idx in list_of_C:
                C_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[C_atom] = new_list
    # Actual Replacement
    COlist.sort()
    rem_atoms = 0
    for number in COlist:
        adj_num = number-rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
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
        adj_num = number - rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        # print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
    return moll, replacements

def repl_atommap_SOOO(moll, replacements:dict):
    '''
    Replace the substructure S(=O)(=O)O by S using RWMol and AtomMapNumbers

    input:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the S that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule

    returns:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the S that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    -- but both are adjusted with the replaced S(=O)(=O)O substructures
    '''
    list_of_S = [Sidx[0] for Sidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('S'))]
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    SOOOlist = []
    query = Chem.MolFromSmiles('S(=O)(=O)O')
    idx_SOOO = moll.GetSubstructMatches(query)
    for idx__SOOO in idx_SOOO:
        new_list = []
        for indiv_idx in idx__SOOO:
            if indiv_idx not in list_of_S:
                if indiv_idx not in SOOOlist:
                    SOOOlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_S:
                S_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[S_atom] = new_list
    # Actual Replacement
    SOOOlist.sort()
    rem_atoms = 0
    for number in SOOOlist:
        adj_num = number - rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        # print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
    return moll, replacements

def repl_atommap_SOO(moll, replacements:dict):
    '''
    Replace the substructure S(=O)(=O) by S using RWMol and AtomMapNumbers

    input:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the S that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule

    returns:
    moll - molecule in MOL format from RDKit also containing AtomMapNumbers
    replacements - dictionary with as key the AtomMapNumber (int) of the S that
    remains in the molecule and as values a list of the AtomMapNumbers (int)
    that are removed from the molecule
    -- but both are adjusted with the replaced S(=O)(=O) substructures
    '''
    list_of_S = [Sidx[0] for Sidx in
                 moll.GetSubstructMatches(Chem.MolFromSmiles('S'))]
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    SOOlist = []
    query = Chem.MolFromSmiles('S(=O)(=O)')
    idx_SOO = moll.GetSubstructMatches(query)
    for idx__SOO in idx_SOO:
        new_list = []
        for indiv_idx in idx__SOO:
            if indiv_idx not in list_of_S:
                if indiv_idx not in SOOlist:
                    SOOlist.append(indiv_idx)
                new_list += [moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()]
            elif indiv_idx in list_of_S:
                S_atom = moll.GetAtomWithIdx(indiv_idx).GetAtomMapNum()
        replacements[S_atom] = new_list
    # Actual Replacement
    SOOlist.sort()
    rem_atoms = 0
    for number in SOOlist:
        adj_num = number - rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        # print(nblist)
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
    return moll, replacements

def repl_atommap_benzene(moll, replacements:dict):
    '''
    Replace the substructure C1CCCCC1 by C using RWMol and AtomMapNumbers

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
    print(Chem.MolToSmiles(moll))
    rwmol = rdkit.Chem.rdchem.RWMol()
    rwmol.InsertMol(moll)
    CClist = []
    query = Chem.MolFromSmiles('C1CCCCC1')
    idx_CC = moll.GetSubstructMatches(query)
    for idx__CC in idx_CC:
        C_atom = idx__CC[0]
        CClist = list(idx__CC[1:])
        replacements[C_atom] = CClist
    # Actual Replacement
    CClist.sort()
    rem_atoms = 0
    for number in CClist:
        adj_num = number-rem_atoms
        repl_atom = moll.GetAtomWithIdx(adj_num)
        nblist = [x.GetAtomMapNum() for x in repl_atom.GetNeighbors()]
        rwmol.RemoveAtom(adj_num)
        if len(nblist) == 1:
            # print(Chem.MolToSmiles(rwmol))
            rem_atoms += 1
            moll = rwmol.GetMol()
            continue
        elif len(nblist) == 2:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            # print(idxlist)
            # print(Chem.MolToSmiles(rwmol))
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        elif len(nblist) == 3:
            idxlist = []
            moll = rwmol.GetMol()
            # print(Chem.MolToSmiles(moll))
            for atom in moll.GetAtoms():
                if atom.GetAtomMapNum() in nblist:
                    # print(atom.GetAtomMapNum(), atom.GetIdx())
                    idxlist.append(atom.GetIdx())
            rwmol.RemoveBond(idxlist[0], idxlist[1])
            rwmol.AddBond(idxlist[0], idxlist[1],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[0], idxlist[2])
            rwmol.AddBond(idxlist[0], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
            rwmol.RemoveBond(idxlist[1], idxlist[2])
            rwmol.AddBond(idxlist[1], idxlist[2],
                          rdkit.Chem.rdchem.BondType.SINGLE)
        else:
            print('STILL GOT A PROBLEM HERE WOOHOO')
            print(len(nblist))
        rem_atoms += 1
        moll = rwmol.GetMol()
        # print(Chem.MolToSmiles(moll))
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
    print(repl_dicts)
    # print(index_dicts)
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

def return_replaced2(repl_dicts:dict, index_dicts:dict):
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
                if value not in subgraph:
                    continue
                index = index_dicts[subgraphlength].index(subgraph)
                for val in repl_dicts[value]:
                    subgraph.add(val)
                index_dicts[subgraphlength][index] = subgraph
    return index_dicts

