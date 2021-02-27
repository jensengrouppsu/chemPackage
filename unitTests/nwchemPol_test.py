import unittest
from chemPackage import collect
import numpy as np
class TestNWExc(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Pol.out')

#    def test_nexcite(self):
#        self.assertEqual(self.test.nexcite, 4, 'incorrect nexcite')
    def test_npol(self):
        self.assertEqual(self.test.nexcite, 1, 'incorrect npol')

#    def test_qmPol(self):
#        pol = np.array([0.255348, 0.336764, 0.336769, 0.45269 ])
#        np.testing.assert_array_almost_equal(self.test.qm_pol, pol)
