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
        3 : ['NNC(', 'NC(O)=', 'NC(N', 'C(O)=N', 'C(NN'], 4 : ['NNC(O)=', 'NC(O)=N', 'NC(NN', 'C(O)=NN'], \
        5: ['NNC(O)=N', 'NC(O)=NN', 'NNC(NN']}
        self.assertEqual(result, expected)    