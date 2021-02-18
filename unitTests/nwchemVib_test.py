import unittest
from chemPackage import collect
import numpy as np
class TestNWVib(unittest.TestCase):
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

    def test_nModes(self):
        self.assertEqual(self.test.nmodes, 12, 'incorrect # of normal modes')
        self.assertEqual(np.shape(self.test.normal_modes), (12,4,3), 'wrong shape for normal modes array')
        nm0 = np.array(((-0.017761,  0.02229 , -0.002359),
                        (0.063445, -0.19627 , -0.312672),
                        (0.200778, -0.058942,  0.307991),
                        (0.013177, -0.008675, -0.002397)))
        nm1 = np.array(((-0.015589, -0.01384 , -0.020262),
                        (0.137969,  0.247559, -0.205675),
                        (-0.180876, -0.071286,  0.008945),
                        (0.254674, -0.187991,  0.135954)))
        nm2 = np.array([[ 0.012879,  0.013271, -0.211822],
                        [-0.061085, -0.052533, -0.086751],
                        [-0.063249, -0.054697, -0.106205],
                        [-0.060293, -0.05549 , -0.44251 ]])
        nm3 = np.array([[-0.01634 , -0.021402, -0.042451],
                        [ 0.033448, -0.024203, -0.154656],
                        [ 0.193481,  0.135824, -0.266708],
                        [-0.02513 ,  0.194401,  0.294005]])
        nm4 = np.array([[-0.002095,  0.255491, -0.001903],
                        [ 0.022595,  0.312314, -0.022914],
                        [-0.042032,  0.247686, -0.008336],
                        [ 0.046243,  0.224038,  0.025551]])
        nm5 = np.array([[ 0.265209, -0.000984, -0.001302],
                        [ 0.264563,  0.006136,  0.004375],
                        [ 0.262693,  0.004277, -0.014459],
                        [ 0.265242,  0.003598,  0.006201]])
        nm6 = np.array([[ 0.      ,  0.      , -0.060425],
                        [-0.113518,  0.030419,  0.279847],
                        [ 0.030419, -0.113518,  0.279847],
                        [ 0.083099,  0.083099,  0.279847]])
        nm7 = np.array([[-0.024831,  0.023958, -0.      ],
                        [-0.11819 , -0.168986,  0.121638],
                        [ 0.181333,  0.120007, -0.11915 ],
                        [ 0.281851, -0.283886, -0.002488]])
        nm8 = np.array([[-0.023958, -0.024831, -0.      ],
                        [ 0.052924,  0.348186,  0.067355],
                        [ 0.341916,  0.048657,  0.071666],
                        [-0.061976, -0.051855, -0.139021]])
        nm9 = np.array([[-0.      , -0.      ,  0.021518],
                        [-0.278711,  0.07468 , -0.099656],
                        [ 0.07468 , -0.278711, -0.099656],
                        [ 0.204031,  0.204031, -0.099656]])
        nm10 = np.array([[ 0.029481, -0.02931 ,  0.      ],
                         [-0.33386 ,  0.08348 , -0.145317],
                         [-0.083121,  0.332761,  0.144832],
                         [ 0.007396, -0.009025,  0.000486]])
        nm11 = np.array([[ 0.02931 ,  0.029481,  0.      ],
                         [-0.187997,  0.0608  , -0.083337],
                         [ 0.061286, -0.189935, -0.084182],
                         [-0.280503, -0.280453,  0.167519]])

        np.testing.assert_array_almost_equal(self.test.normal_modes[0,:,:], nm0)
        np.testing.assert_array_almost_equal(self.test.normal_modes[1,:,:], nm1)
        np.testing.assert_array_almost_equal(self.test.normal_modes[2,:,:], nm2)
        np.testing.assert_array_almost_equal(self.test.normal_modes[3,:,:], nm3)
        np.testing.assert_array_almost_equal(self.test.normal_modes[4,:,:], nm4)
        np.testing.assert_array_almost_equal(self.test.normal_modes[5,:,:], nm5)
        np.testing.assert_array_almost_equal(self.test.normal_modes[6,:,:], nm6)
        np.testing.assert_array_almost_equal(self.test.normal_modes[7,:,:], nm7)
        np.testing.assert_array_almost_equal(self.test.normal_modes[8,:,:], nm8)
        np.testing.assert_array_almost_equal(self.test.normal_modes[9,:,:], nm9)
        np.testing.assert_array_almost_equal(self.test.normal_modes[10,:,:], nm10)
        np.testing.assert_array_almost_equal(self.test.normal_modes[11,:,:], nm11)
    
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