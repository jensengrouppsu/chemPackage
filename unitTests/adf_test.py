import unittest
from chemPackage import collect

class TestADFMethods(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/raman-benzene.out')

    def test_init(self):
        self.assertEqual(type(self.test.program), type(str('ADF')), 'incorrect program')
        self.assertEqual(self.test.program, 'ADF', 'incorrect program')
        self.assertEqual(self.test.project, 'all', 'incorrect project')
        self.assertEqual(self.test.filetype, 'out', 'incorrect filetype')
        self.assertEqual(self.test.filename, 'unitTests/raman-benzene.out', 'incorrect name')
    def test_calctype(self):
        ct_array = {'DFT', 'FREQUENCIES', 'POLARIZABILITY', 'RAMAN', 'STATIC'}
        self.assertEqual(type(self.test.calctype), type(ct_array), 'incorrect calctype type')
        self.assertEqual(len(self.test.calctype), len(ct_array), 'mismatch length')
        self.assertEqual(self.test.calctype, ct_array, 'mismatch string')
    def test_collect(self):
        # unit tests for testing things that must be correct in adf collection
        self.assertEqual(self.test._abort, True, 'incorrect abort setting')
if __name__ == '__main__':
    unittest.main()