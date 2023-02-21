#!/usr/bin/env python3
import unittest 

from GraphMiner import increment
from GraphMiner import select_on_size, list_nodes, rdkit_parse, breadth_fs, subgraphs_smiles
    # def test_increment(self):
    #     result = increment(9)
    #     expected = 10
    #     self.assertEqual(result, expected)

class TestSequence(unittest.TestCase):   
    def test_subgraph_mining(self):
        mol_list = ['C', 'N', 'O', 'Br', 'S']
        selected_smile = select_on_size('NNC(O)=NN')
        list_node = list_nodes(selected_smile, mol_list)
        dictnode = rdkit_parse(selected_smile, list_node)
        subgraphdict = breadth_fs(list_node, dictnode)
        result = subgraphs_smiles(subgraphdict, selected_smile)
        expected = {1 : ['N', 'N', 'C(', 'O)=', 'N', 'N'], 2 : ['NN', 'NC(', 'C(O)=', 'C(N', 'NN'], \
        3 : ['NNC(', 'NC(O)=', 'NC(N', 'C(O)=N', 'C(NN'], 4 : ['NNC(O)=', 'NNC(N', 'NC(O)=N', 'NC(NN', 'C(O)=NN'], \
        5: ['NNC(O)=N', 'NNC(NN', 'NC(O)=NN'], 6 : ['NNC(O)=NN']}
        self.assertEqual(result, expected)

    def test_partly_subgraphs(self):
        mol_list = ['C', 'N', 'O', 'Br', 'S']
        selected_smile = select_on_size('NNC(O)=NN')
        list_node = list_nodes(selected_smile, mol_list)
        dictnode = rdkit_parse(selected_smile, list_node)
        subgraphdict = breadth_fs(list_node, dictnode)
        result = subgraphs_smiles(subgraphdict, selected_smile)
        expected = ['NNC(O)=', 'NNC(N', 'NC(O)=N', 'NC(NN', 'C(O)=NN']
        self.assertIn(expected, result.values())
        #Checks if expected is in result

    def test_ring_structure(self):
        mol_list = ['C', 'N', 'O', 'Br', 'S']
        selected_smile = select_on_size('CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=C(C=C3)O)N)C(=O)O)C')
        list_node = list_nodes(selected_smile, mol_list)
        dictnode = rdkit_parse(selected_smile, list_node)
        subgraphdict = breadth_fs(list_node, dictnode)
        result = subgraphs_smiles(subgraphdict, selected_smile)
        expected = 'C1(C(N2C(S1)C(C2=O)'
        self.assertIn(expected, result[8])
