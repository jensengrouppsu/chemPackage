import unittest
from chemPackage import collect

class TestNWClass(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Geom.out')
        self.test2 = collect('unitTests/nh3Freq.out')
        self.test3 = collect('unitTests/nh3Pol.out')
        self.test4 = collect('unitTests/nh3Pol2.out')

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
        
    def test_dipole(self):
        self.assertEqual(self.test.dipole[0], -1.6998385895583998, 'incorrect dipole x')
        self.assertEqual(self.test.dipole[1], 0.44072101619589993, 'incorrect dipole y')
        self.assertEqual(self.test.dipole[2], -0.0, 'incorrect dipole z')

    def test_symmetry(self):
        self.assertEqual(self.test.symmetry, 'Cs', 'incorrect symmetry')

    def test_calctype(self):
        self.assertEqual(self.test.calctype, {'RESTRICTED', 'DFT', 'GEOMETRY'}, 'incorrect symmetry')
        self.assertEqual(self.test2.calctype, {'RESTRICTED', 'DFT', 'FREQUENCIES'}, 'incorrect symmetry')
        self.assertEqual(self.test3.calctype, {'RESTRICTED', 'DFT', 'OPTICAL ROTATION', 'POLARIZABILITY', 'FD'}, 'incorrect symmetry')
        self.assertEqual(self.test4.calctype, {'RESTRICTED', 'DFT', 'OPTICAL ROTATION', 'POLARIZABILITY', 'STATIC'}, 'incorrect symmetry')