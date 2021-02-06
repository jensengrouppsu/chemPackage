import unittest
from chemPackage import collect
import numpy as np

class TestADFMethods(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/raman-benzene.out')

    def test_vFrequencies(self):
        self.assertEqual(len(self.test.v_frequencies), 30, 'incorrect vFreq length')
        self.assertEqual(self.test.v_frequencies[0], 411.02, 'incorrect low freq')
        self.assertEqual(self.test.v_frequencies[17], 1134.248, 'incorrect mid freq')
        self.assertEqual(self.test.v_frequencies[29], 3137.143, 'incorrect high freq')

if __name__ == '__main__':
    unittest.main()