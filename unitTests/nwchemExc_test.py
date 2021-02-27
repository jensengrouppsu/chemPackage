import unittest
from chemPackage import collect
import numpy as np
class TestNWExc(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Exc.out')

    def test_nexcite(self):
        self.assertEqual(self.test.nexcite, 4, 'incorrect nexcite')

    def test_ExcE(self):
        ExcE = np.array([0.255348, 0.336764, 0.336769, 0.45269 ])
        np.testing.assert_array_almost_equal(self.test.excitation_energies, ExcE)

    def test_osc(self):
        Osc = np.array([0.037808, 0.031868, 0.031864, 0.19002 ])
        np.testing.assert_array_almost_equal(self.test.oscillator_strengths, Osc)

    def test_ExcSym(self):
        Sym = np.array(['A', 'A', 'A', 'A'])
        np.testing.assert_array_equal(self.test.excitation_symmetries, Sym)

    def test_ExcType(self):
        Sym = np.array(['SS', 'SS', 'SS', 'SS'])
        np.testing.assert_array_equal(self.test.excitation_type, Sym)

    def test_TDM(self):
        TDM = np.array([[ 4.5619e-01, -1.1826e-01, -0.0000e+00],
                        [-9.4550e-02, -3.6470e-01, -0.0000e+00],
                        [ 0.0000e+00,  0.0000e+00, -3.7673e-01],
                        [ 0.0000e+00,  2.0000e-05, -7.9350e-01]])
        np.testing.assert_array_almost_equal(self.test.TDM, TDM)

    def test_TQM(self):
        TQM = np.array([[-2.5638e-01, -1.4990e-02,  0.0000e+00, -2.5638e-01, -1.4990e-02, 0.0000e+00],
                        [ 2.3376e-01,  4.4838e-01,  0.0000e+00,  2.3376e-01,  4.4838e-01, 0.0000e+00],
                        [-0.0000e+00, -0.0000e+00,  3.8467e-01, -0.0000e+00, -0.0000e+00, 3.8467e-01],
                        [-1.0000e-05, -1.0000e-05,  3.2939e-01, -1.0000e-05, -1.0000e-05, 3.2939e-01]])
        np.testing.assert_array_almost_equal(self.test.TQM, TQM)