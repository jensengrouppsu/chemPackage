import unittest
from chemPackage import collect
import numpy as np
from numpy.testing import assert_array_almost_equal

class TestADFMethods(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_pol.out')

    def test_eFrequencies(self):
        self.assertEqual(np.shape(self.test.e_frequencies), (1,), 'wrong shape for e frequencies')
    def test_nPol(self):
        self.assertEqual(self.test.npol, 1, 'wrong number of polarizations')

    def test_qmPol(self):
        pol = np.array([[[7.80574+0.40954j, -0.+0.j   , 0. + 0.j],
                        [-0.+0.j   ,9.04814+0.17714j, 0. +0.j],
                        [-0.+0.j, -0.+0.j, 8.31199+0.25378j]]])
        np.testing.assert_array_almost_equal(self.test.qm_pol, pol)