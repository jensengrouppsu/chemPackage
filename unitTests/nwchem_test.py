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

    def test_energy(self):
        self.assertEqual(self.test.energy['total DFT'], -56.509475409542, 'incorrect energy')
        self.assertEqual(self.test.energy['one-electron'], -99.473135904705, 'incorrect energy')
        self.assertEqual(self.test.energy['Coulomb'], 39.146396077606, 'incorrect energy')
        self.assertEqual(self.test.energy['XC'], -8.026930069048, 'incorrect energy')
        self.assertEqual(self.test.energy['nuc. repulsion'], 11.844194486605, 'incorrect energy')
        self.assertEqual(self.test.energy['total'], -56.509475409542, 'incorrect energy')
