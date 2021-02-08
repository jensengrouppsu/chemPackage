import unittest
from chemPackage import collect
import numpy as np
from numpy.testing import assert_array_almost_equal

class TestADFMethods(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_freq.out')