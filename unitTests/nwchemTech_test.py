import unittest
from chemPackage import collect
import numpy as np
from datetime import datetime, timedelta
class TestNWClass(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Freq.out')

    def test_nproc(self):
        self.assertEqual(self.test.nprocs, 1, 'incorrect nproc')

    def test_startTime(self):
        self.assertEqual(self.test.start, datetime(2021, 2, 14, 17, 41, 14), 'incorrect start time')
        
    def test_termination(self):
        self.assertEqual(self.test.termination, 'NORMAL TERMINATION', 'incorrect termination')

    def test_host(self):
        self.assertEqual(self.test.host, 'neon.science.psu.edu', 'wrong host')