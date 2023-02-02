#!/usr/bin/env python3
import unittest 

from GraphMiner import increment


class TestSequence(unittest.TestCase):

    def test_increment(self):
        result = increment(9)
        expected = 10
        self.assertEqual(result, expected)