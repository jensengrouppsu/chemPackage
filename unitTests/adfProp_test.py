import unittest
from chemPackage import collect
import numpy as np

class TestADFMethods(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_freq.out')

    def test_vFrequencies(self):
        self.assertEqual(len(self.test.v_frequencies), 3, 'incorrect vFreq length')
        self.assertEqual(self.test.v_frequencies[0], 1578.216, 'incorrect low freq')
        self.assertEqual(self.test.v_frequencies[1], 3616.061, 'incorrect mid freq')
        self.assertEqual(self.test.v_frequencies[2], 3727.173, 'incorrect high freq')

if __name__ == '__main__':
    unittest.main()