import unittest
from chemPackage import collect
import numpy as np
class TestNWExc(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Pol2.out')

#    def test_nexcite(self):
#        self.assertEqual(self.test.nexcite, 4, 'incorrect nexcite')
    def test_npol(self):
        self.assertEqual(self.test.npol, 1, 'incorrect npol')

    def test_qmPol(self):
        pol = np.array([[[ 7.357008,  0.727041, -0.      ],
                         [ 0.727041,  9.971928,  0.      ],
                         [-0.      ,  0.      , 10.159759]]])
        np.testing.assert_array_almost_equal(self.test.qm_pol, pol)
