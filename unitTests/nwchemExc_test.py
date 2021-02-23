import unittest
from chemPackage import collect
import numpy as np
class TestNWExc(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Exc.out')

    def test_nexcite(self):
        self.assertEqual(self.test.nexcite, 4, 'incorrect nexcite')