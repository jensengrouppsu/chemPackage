import unittest
from chemPackage import collect

class TestADFMethods(unittest.TestCase):
    def setUp(self):
        self.test = collect('../adf/raman/raman-benzene.orig.out')

    def test_calctype(self):
        ct_array = {'DFT', 'FREQUENCIES', 'POLARIZABILITY', 'RAMAN', 'STATIC'}
        self.assertEqual(type(self.test.calctype), type(ct_array), 'incorrect calctype type')
        self.assertEqual(len(self.test.calctype), len(ct_array), 'mismatch length')
        self.assertEqual(self.test.calctype, ct_array, 'mismatch string')
        

if __name__ == '__main__':
    unittest.main()