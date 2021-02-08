import unittest
from chemPackage import collect
import numpy as np
from numpy.testing import assert_array_almost_equal

class TestADFMethods(unittest.TestCase):
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