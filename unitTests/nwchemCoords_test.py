import unittest
from chemPackage import collect
import numpy as np
class TestNWCoords(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Geom.out')

    def test_natoms(self):
        self.assertEqual(self.test.natoms, 4, 'incorrect natoms')

    def test_coordinates(self):
        coords = np.array([[ 0.219725,  0.07246 ,  0.      ],
                           [-0.057493,  0.628731, -0.812078],
                           [-0.057493,  0.628731,  0.812078],
                           [-0.410566, -0.73287 ,  0.      ]])
        np.testing.assert_array_almost_equal(self.test.coordinates, coords)