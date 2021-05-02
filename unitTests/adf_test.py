import unittest
from chemPackage import collect
import numpy as np
from numpy.testing import assert_array_almost_equal

class TestADFMethods(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/raman-benzene.out')

    def test_init(self):
        self.assertEqual(type(self.test.program), type(str('ADF')), 'incorrect program')
        self.assertEqual(self.test.program, 'ADF', 'incorrect program')
        self.assertEqual(self.test.project, 'all', 'incorrect project')
        self.assertEqual(self.test.filetype, 'out', 'incorrect filetype')
        self.assertEqual(self.test.filename, 'unitTests/raman-benzene.out', 'incorrect name')
    def test_calctype(self):
        ct_array = {'DFT', 'FREQUENCIES', 'POLARIZABILITY', 'RAMAN', 'STATIC'}
        self.assertEqual(type(self.test.calctype), type(ct_array), 'incorrect calctype type')
        self.assertEqual(len(self.test.calctype), len(ct_array), 'mismatch length')
        self.assertEqual(self.test.calctype, ct_array, 'mismatch string')
    def test_collect(self):
        # unit tests for testing things that must be correct in adf collection
        self.assertEqual(self.test._abort, True, 'incorrect abort setting')

class TestADFOrbitals(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_pol.out')

    def test_Orbital(self):
        orbE = np.array([-25.34, -13.259, -9.557, -7.398, -0.228, 2.232])
        np.testing.assert_array_almost_equal(self.test.orbital_energies, orbE)
        occ = np.array([2., 2., 2., 2., 0., 0.])
        np.testing.assert_array_almost_equal(self.test.orbital_occ, occ)
        spin = np.array([0, 0, 0, 0, 0, 0])
        np.testing.assert_array_almost_equal(self.test.orbital_spin, spin)

class TestADFPol(unittest.TestCase):
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

class TestADFVib(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_freq.out')

    def test_vFrequencies(self):
        self.assertEqual(len(self.test.v_frequencies), 3, 'incorrect vFreq length')
        self.assertEqual(self.test.v_frequencies[0], 1578.216, 'incorrect low freq')
        self.assertEqual(self.test.v_frequencies[1], 3616.061, 'incorrect mid freq')
        self.assertEqual(self.test.v_frequencies[2], 3727.173, 'incorrect high freq')

    def test_nModes(self):
        self.assertEqual(self.test.nmodes, 3, 'incorrect # of normal modes')
        nm1 = np.array((((0. , -0.22394376, -0.29859169), 
                       (0. ,  0.22394376, -0.29859169),  
                       (0. ,  0. ,          0.03758867)),
                       ((0., -0.31096703, 0.20731136),
                        (0.,  0.31096703, 0.20731136),
                        (0.,  0.,        -0.02591392)),
                       ((0.,  0.2935258, -0.23058963),
                        (0.,  0.2935258,  0.23058963),
                        (0., -0.03702127, 0.))))
        self.assertEqual(np.shape(self.test.normal_modes), (3,3,3), 'wrong shape for normal modes array')
        np.testing.assert_array_almost_equal(self.test.normal_modes, nm1)
    
    def test_IR(self):
        self.assertEqual(len(self.test.IR), 3, 'incorrect # IR absorption intensities')
        self.assertEqual(self.test.IR[0], 80.448586, 'incorrect first IR intensity')
        self.assertEqual(self.test.IR[1], 6.537245, 'incorrect second IR intensity')
        self.assertEqual(self.test.IR[2], 89.340161, 'incorrect third IR intensity')

if __name__ == '__main__':
    unittest.main()