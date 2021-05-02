import unittest
from chemPackage import collect
import numpy as np
from datetime import datetime, timedelta
from numpy.testing import assert_array_almost_equal
from chemPackage.coords import Coordinates

class TestNWClass(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Geom.out')
        self.test2 = collect('unitTests/nh3Freq.out')
        self.test3 = collect('unitTests/nh3Pol.out')
        self.test4 = collect('unitTests/nh3Pol2.out')
        self.test5 = collect('unitTests/DIM_spe.out')

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
        self.assertEqual(self.test5.calctype, {'RESTRICTED', 'DFT', 'DIM'}, 'incorrect symmetry')

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

    def test_atoms(self):
        atoms = np.array(['N', 'H', 'H', 'H'], dtype='<U1')
        np.testing.assert_array_equal(self.test.atoms, atoms)

    def test_elements(self):
        ele = {'N', 'H'}
        self.assertEqual(self.test.elements, ele)

    def test_nelements(self):
        self.assertEqual(self.test.nelements, 2, 'wrong length of elements')

    def test_optGeom(self):
        initCoords = np.array([[ 0.050971, -0.099509,  0.      ],
                               [-0.055838,  0.599485, -0.5     ],
                               [-0.055838,  0.599485,  0.5     ],
                               [-0.245121, -0.50241 ,  0.      ]], dtype=np.float64)
        np.testing.assert_array_almost_equal(self.test.initial_geometry, initCoords)

    def test_gsGradient(self):
        grad = np.array([[-1.7e-05,  3.1e-05,  0.0e+00],
                         [ 1.2e-05,  2.0e-06,  0.0e+00],
                         [ 1.2e-05,  2.0e-06, -0.0e+00],
                         [-7.0e-06, -3.5e-05,  0.0e+00]], dtype=np.float64)
        np.testing.assert_array_almost_equal(self.test.gs_gradient['total DFT gradient'], grad)

class TestNWTech(unittest.TestCase):
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

class TestNWExc(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Exc.out')

    def test_nexcite(self):
        self.assertEqual(self.test.nexcite, 4, 'incorrect nexcite')

    def test_ExcE(self):
        ExcE = np.array([0.255348, 0.336764, 0.336769, 0.45269 ])
        np.testing.assert_array_almost_equal(self.test.excitation_energies, ExcE)

    def test_osc(self):
        Osc = np.array([0.037808, 0.031868, 0.031864, 0.19002 ])
        np.testing.assert_array_almost_equal(self.test.oscillator_strengths, Osc)

    def test_ExcSym(self):
        Sym = np.array(['A', 'A', 'A', 'A'])
        np.testing.assert_array_equal(self.test.excitation_symmetries, Sym)

    def test_ExcType(self):
        Sym = np.array(['SS', 'SS', 'SS', 'SS'])
        np.testing.assert_array_equal(self.test.excitation_type, Sym)

    def test_TDM(self):
        TDM = np.array([[ 4.5619e-01, -1.1826e-01, -0.0000e+00],
                        [-9.4550e-02, -3.6470e-01, -0.0000e+00],
                        [ 0.0000e+00,  0.0000e+00, -3.7673e-01],
                        [ 0.0000e+00,  2.0000e-05, -7.9350e-01]])
        np.testing.assert_array_almost_equal(self.test.TDM, TDM)

    def test_TQM(self):
        TQM = np.array([[-2.5638e-01, -1.4990e-02,  0.0000e+00, -2.5638e-01, -1.4990e-02, 0.0000e+00],
                        [ 2.3376e-01,  4.4838e-01,  0.0000e+00,  2.3376e-01,  4.4838e-01, 0.0000e+00],
                        [-0.0000e+00, -0.0000e+00,  3.8467e-01, -0.0000e+00, -0.0000e+00, 3.8467e-01],
                        [-1.0000e-05, -1.0000e-05,  3.2939e-01, -1.0000e-05, -1.0000e-05, 3.2939e-01]])
        np.testing.assert_array_almost_equal(self.test.TQM, TQM)

class TestNWPol(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/nh3Pol2.out')
        self.testDamp = collect('unitTests/nh3Pol.out')

#    def test_nexcite(self):
#        self.assertEqual(self.test.nexcite, 4, 'incorrect nexcite')
    def test_npol(self):
        self.assertEqual(self.test.npol, 1, 'incorrect npol')
        self.assertEqual(self.testDamp.npol, 1, 'incorrect damping npol')

    def test_qmPol(self):
        pol = np.array([[[ 7.357008,  0.727041, -0.      ],
                         [ 0.727041,  9.971928,  0.      ],
                         [-0.      ,  0.      , 10.159759]]])
        np.testing.assert_array_almost_equal(self.test.qm_pol, pol)

        polDamp = np.array([[[ 7.356923+0.010561j,  0.727051-0.000624j, -0.      +0.j      ],
                             [ 0.727051-0.000624j,  9.971879+0.008315j,  0.      +0.j      ],
                             [-0.      +0.j      ,  0.      +0.j      , 10.159713+0.008153j]]])
        np.testing.assert_array_almost_equal(self.testDamp.qm_pol, polDamp)

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

class TestADFOrbitals(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_pol.out')

    def test_Orbital(self):
        orbE = np.array([-25.34, -13.259, -9.557, -7.398, -0.228, 2.232])
        np.testing.assert_array_almost_equal(self.test.orbital_energies, orbE)
        occ = np.array([2., 2., 2., 2., 0., 0.])
        np.testing.assert_array_almost_equal(self.test.orbital_occ, occ)
        spin = np.array([0, 0, 0, 0, 0, 0])
        np.testing.assert_array_almost_equal(self.test.orbital_spin, spin)

class TestADFPol(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_pol.out')

    def test_eFrequencies(self):
        self.assertEqual(np.shape(self.test.e_frequencies), (1,), 'wrong shape for e frequencies')
    def test_nPol(self):
        self.assertEqual(self.test.npol, 1, 'wrong number of polarizations')

    def test_qmPol(self):
        pol = np.array([[[7.80574+0.40954j, -0.+0.j   , 0. + 0.j],
                        [-0.+0.j   ,9.04814+0.17714j, 0. +0.j],
                        [-0.+0.j, -0.+0.j, 8.31199+0.25378j]]])
        np.testing.assert_array_almost_equal(self.test.qm_pol, pol)

class TestADFVib(unittest.TestCase):
    def setUp(self):
        self.test = collect('unitTests/h2o_freq.out')

    def test_vFrequencies(self):
        self.assertEqual(len(self.test.v_frequencies), 3, 'incorrect vFreq length')
        self.assertEqual(self.test.v_frequencies[0], 1578.216, 'incorrect low freq')
        self.assertEqual(self.test.v_frequencies[1], 3616.061, 'incorrect mid freq')
        self.assertEqual(self.test.v_frequencies[2], 3727.173, 'incorrect high freq')

    def test_nModes(self):
        self.assertEqual(self.test.nmodes, 3, 'incorrect # of normal modes')
        nm1 = np.array((((0. , -0.22394376, -0.29859169), 
                       (0. ,  0.22394376, -0.29859169),  
                       (0. ,  0. ,          0.03758867)),
                       ((0., -0.31096703, 0.20731136),
                        (0.,  0.31096703, 0.20731136),
                        (0.,  0.,        -0.02591392)),
                       ((0.,  0.2935258, -0.23058963),
                        (0.,  0.2935258,  0.23058963),
                        (0., -0.03702127, 0.))))
        self.assertEqual(np.shape(self.test.normal_modes), (3,3,3), 'wrong shape for normal modes array')
        np.testing.assert_array_almost_equal(self.test.normal_modes, nm1)
    
    def test_IR(self):
        self.assertEqual(len(self.test.IR), 3, 'incorrect # IR absorption intensities')
        self.assertEqual(self.test.IR[0], 80.448586, 'incorrect first IR intensity')
        self.assertEqual(self.test.IR[1], 6.537245, 'incorrect second IR intensity')
        self.assertEqual(self.test.IR[2], 89.340161, 'incorrect third IR intensity')

class TestCoordFuncs(unittest.TestCase):

    def test_dtoLine(self):
        f = Coordinates()
        a = np.array([0, 0, 0])
        b = np.array([0, 0, 5])
        c = np.array([2, 0, 0])
        answer = f.dtoLine(a, b, c)

        self.assertEqual(answer, 2)

    def test_cutCylinderCoords(self):
        f = collect('unitTests/sampleCoords.xyz')
        index, atomType, cylinder = f.cutCylinder(4, 5, radius = 5.0)
        answer = np.array([[0.0, 0.0, 0.0],
                                 [0.75, 0.75, 0.0],
                                 [-0.75, 0.75, 0.0],
                                 [-5.0, 0.0, 0.0],
                                 [5.0, 0.0, 0.0]
                                               ])
        np.testing.assert_array_almost_equal(cylinder, answer)

    def test_cutCylinderIndex(self):
        f = collect('unitTests/sampleCoords.xyz')
        index, atomType, cylinder = f.cutCylinder(4, 5, radius = 5.0)
        answer = np.array([1, 2, 3, 4, 5])
        np.testing.assert_array_almost_equal(index, answer)

    def test_cutCylinderAtomType(self):
        f = collect('unitTests/sampleCoords.xyz')
        index, atomType, cylinder = f.cutCylinder(4, 5, radius = 5.0)
        answer = np.array(['H', 'O', 'O', 'Ag', 'Ag'])
        np.testing.assert_array_equal(atomType, answer)

if __name__ == '__main__':
    unittest.main()