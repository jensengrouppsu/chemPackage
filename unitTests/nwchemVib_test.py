import unittest
from chemPackage import collect
import numpy as np
class TestNWClass(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Freq.out')

    def test_vFrequencies(self):
        self.assertEqual(len(self.test.v_frequencies), 12, 'incorrect vFreq length')
        self.assertEqual(self.test.v_frequencies[0], -0.0, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[1], -0.0, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[2], -0.0, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[3], -0.0, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[4], -0.0, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[5], -0.0, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[6], 1061.43, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[7], 1656.51, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[8], 1656.51, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[9], 3444.4, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[10], 3561.84, 'incorrect freq')
        self.assertEqual(self.test.v_frequencies[11], 3561.84, 'incorrect freq')

#    def test_nModes(self):
#        self.assertEqual(self.test.nmodes, 12, 'incorrect # of normal modes')
#        nm1 = np.array((((0. , -0.22394376, -0.29859169), 
#                       (0. ,  0.22394376, -0.29859169),  
#                       (0. ,  0. ,          0.03758867)),
#                       ((0., -0.31096703, 0.20731136),
#                        (0.,  0.31096703, 0.20731136),
#                        (0.,  0.,        -0.02591392)),
#                       ((0.,  0.2935258, -0.23058963),
#                        (0.,  0.2935258,  0.23058963),
#                        (0., -0.03702127, 0.))))
#        self.assertEqual(np.shape(self.test.normal_modes), (12,4,3), 'wrong shape for normal modes array')
#        np.testing.assert_array_almost_equal(self.test.normal_modes, nm1)
    
    def test_IR(self):
        self.assertEqual(len(self.test.IR), 12, 'incorrect # IR absorption intensities')
        self.assertEqual(self.test.IR[0], 64.974, 'incorrect first IR intensity')
        self.assertEqual(self.test.IR[1], 20.092, 'incorrect second IR intensity')
        self.assertEqual(self.test.IR[2], 9.088, 'incorrect third IR intensity')
        self.assertEqual(self.test.IR[3], 55.21, 'incorrect IR intensity(4)')
        self.assertEqual(self.test.IR[4], 0.108, 'incorrect IR intensity(5)')
        self.assertEqual(self.test.IR[5], 0.022, 'incorrect IR intensity(6)')
        self.assertEqual(self.test.IR[6], 118.879, 'incorrect IR intensity(7)')
        self.assertEqual(self.test.IR[7], 8.809, 'incorrect IR intensity(8)')
        self.assertEqual(self.test.IR[8], 8.809, 'incorrect IR intensity(9)')
        self.assertEqual(self.test.IR[9], 0.321, 'incorrect IR intensity(10)')
        self.assertEqual(self.test.IR[10], 0.002, 'incorrect IR intensity(11)')

if __name__ == '__main__':
    unittest.main()