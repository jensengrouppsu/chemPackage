import unittest
from chemPackage import collect

class TestNWClass(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Geom.out')

    def test_init(self):
        self.assertEqual(type(self.test.program), type(str('NWChem')), 'incorrect program')
        self.assertEqual(self.test.program, 'NWChem', 'incorrect program')
        self.assertEqual(self.test.filetype, 'out', 'incorrect filetype')
        self.assertEqual(self.test.filename, 'unitTests/nh3Geom.out', 'incorrect name')