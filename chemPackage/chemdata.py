######################################################################
# Watch out, most of the coding in here is a little more low-level
# than the rest of the chem package, which is a necessary evil to give
# the functionality that we want.
######################################################################

from __future__ import print_function, division
from .errorclass import ChemDataError
from .coords import Coordinates
from .pol import Polarizability
from .excite_mo import Excitations_MOs
from .raman_IR import Raman_IR
from .general import General
from .input import InputFiles
from .charges import charges

# Create the ChemData class to hold and manipulate data
# It is a 'new-style" class.
class ChemData(Coordinates, Polarizability, Excitations_MOs,
               Raman_IR, General, InputFiles):

    def __init__(self, name=None):
        '''Initiallize the ChemData class. This class implements basic
        collection for XYZ files which should be replaced with a
        superclass's own collection implementation.

        '''

        # Initiallize each attribute to none.
        # The __dict__ notation is necessary here bypass the __setattr__ method
        for attr in attrset | privates:
            self.__dict__[attr] = None

        # Initiallize a few to something other than None
        self.__dict__['_privates'] = privates
        self.__dict__['_attrset'] = attrset
        self._attrlist = attrlist
        self._managed_methods = managed_methods
        self._computer_names = computer_names
        self.calctype = set()
        self.subkey = set()
        self.key = {}
        self._i = 0
        self.charges = charges

        # Save the .xyz file if given.
        if name is not None:
            from os.path import splitext
            ftype = splitext(name)[1]
            if ftype != '.xyz':
                raise ValueError (ftype+' not a valid XYZ extention')
            self.filetype = ftype[1:]
            self.filename = name

    def _collect(self, abort=False):
        '''Collect the XYZ data from file.'''

        with open(self.filename) as fl:
            from numpy import array
            # 1st line is # of atoms
            self.natoms = int(fl.readline().strip())
            # 2nd line is title, or is blank
            second = fl.readline().strip()
            self.title = second if second else None
            # Collect remaining lines in an array
            atm_coors = array([ln.strip().split() for ln in fl])
            # Atoms are first index, coordinates are remaining indices
            self.atoms = atm_coors[:,0]
            self.coordinates = array(atm_coors[:,1:], dtype=float)
            self.elements = set(self.atoms)
            self.nelements = len(self.elements)

        # Make sure that the number of atoms and natoms matches
        if self.natoms != len(self.atoms):
            str1 = 'natoms = {0}, len(atoms) = {1}'
            str2 = 'Number of atoms does not match integer at head of file'
            self._raise_or_pass(str1.format(self.natoms, len(self.atoms)),
                                str2, abort)

    def __str__(self):
        '''Defines what is printed when the class is printed.'''
        return 'ChemData object for program {0} and file {1}'.format(
                                                                 self.program,
                                                                 self.filename)

    def __repr__(self):
        '''Defines what is given when the class is given on command-line.'''
        return self.__str__() + ', {0} data fields filled'.format(len(self))

    def __setattr__(self, attr, value):
        '''Prevent creating attributes not defined here.'''
        if attr in self._attrset | self._privates:
            self.__dict__[attr] = value
        else:
            raise AttributeError('Addition of attributes that are not pre-'
                                   'defined is not allowed! (' + attr + ')')

    def __len__(self):
        '''The length of this class is the number of fields filled.'''
        return len([x for x in self])

    def __iter__(self):
        '''Initiallizes iteration'''
        self._i = 0
        return self

    def __contains__(self, key):
        '''Boolean for if a field is filled or not.'''
        if key in ('key', 'subkey', 'calctype'):
            if getattr(self, key):
                return True
            else:
                return False
        elif key not in self._attrset:
            return False
        elif getattr(self, key) is None:
            return False
        else:
            return True

    def next(self):
        '''Looping over the class returns non-empty keys.'''
        # Stop when no more attributes.
        if self._i == len(self._attrlist):
            raise StopIteration
        # Check if the key is None.  If so, find one that isn't.
        key = self._attrlist[self._i]
        while key not in self:
            self._i += 1
            if self._i == len(self._attrlist):
                raise StopIteration
            key = self._attrlist[self._i]
        self._i += 1
        # Return first key found that isn't False.
        return key

    def copy(self):
        '''Returns a copy of the current instance.
        Equivalent to deepcopy(self).

        '''
        from copy import deepcopy
        return deepcopy(self)

    def _raise_or_pass(self, msg):
        '''Subroutine to pass or raise on collection error.
        Makes termination more specific by adding error message.

        '''
        if self.termination is None or 'ERROR:' not in self.termination:
            self.termination = msg
        if self._abort:
            from .errorclass import CollectionError
            raise CollectionError(msg)
        else:
            pass

    def filled(self):
        '''Returns a set of all the attributes on the attribute list that
        are not empty.

        '''
        return set([x for x in self])

    def dump(self):
        '''Method to dump all the collected data to screen.'''

        # Loop over the attribute list and print off what appears in the list
        for k in self:
            print(k + ': ', getattr(self, k), sep='\n')
        for k in self._managed_methods:
            try:
                data = getattr(self, k)
            except (TypeError, AssertionError):
                continue
            if data is None: continue
            print(k + ': ', data, sep='\n')

##################
# Class attributes
##################

attrlist = (
    'program',
    'filetype',
    'filename',
    'host',
    'nprocs',
    'start',
    'real_time',
    'cpu_time',
    'routine_times',
    'termination',
    'title',
    'calctype',
    'key',
    'subkey',
    'mixedkeys',
    'blockkeys',
    'linekeys',
    'singlekeys',
    'natoms',
    'atoms',
    'elements',
    'nelements',
    'coordinates',
    'initial_geometry',
    'symmetry',
    'energy',
    'charges',
    'qmcharge',
    'dipole',
    'nmos',
    'nocc',
    'nvir',
    'orbital_energies',
    'orbital_ids',
    'orbital_occ',
    'orbital_spin',
    'atomic_orbitals',
    'nexcite',
    'excitation_energies',
    'oscillator_strengths',
    'excitation_symmetries',
    'excitation_type',
    'TDM',
    'excTDM',
    'MDM',
    'TQM',
    'transitions',
    'STPM',
    'linear_tpa_strengths',
    'linear_sigma_tpa',
    'circular_tpa_strengths',
    'circular_sigma_tpa',
    'tpa_polarization_ratio',
    'T3PM',
    'linear_3pa_strengths',
    'linear_sigma_3pa',
    'circular_3pa_strengths',
    'circular_sigma_3pa',
    '3pa_polarization_ratio',
    'gs_gradient',
    'es',
    'es_dipole',
    'es_gradient',
    'deltas',
    'dgdip',
    'dtdip',
    'HOMO',
    'LUMO',
    'SOMO',
    'nmodes',
    'normal_modes',
    'IR',
    'v_frequencies',
    'npol',
    'nhpol',
    'e_frequencies',
    'b_e_frequencies',
    'c_e_frequencies',
    'd_e_frequencies',
    'ord',
    'qm_pol',
    'polarizability_atomic',
    'polarizability_derivatives',
    'alpha_derivatives',
    'hyperpolarizability',
    'secondhyperpolarizability',
    'dhpol',
    'dshpol',
    'hyperinvarients',
    'ndim',
    'ndim_elements',
    'dim_atoms',
    'dim_elements',
    'dim_coordinates',
    'dim_cnradii',
    'dim_charges',
    'dim_dipoles',
    'dim_dipole_tot',
    'dim_pol',
    'dim_efficiencies',
    'dim_cross_sections',
    'dim_energy',
    'dimqm_energy',
    'dipole',
    'charges',
    'orbitalplot',
    'orbitaldim',
    'dipfitdens',
    'quadrupole',
    'quadfitdens',
    'atensor',
    'astensor',
    'atensor_derivatives',
    'vroa_intensities',
    'spin_multiplicity',
    'gtensor',
    'gstensor',
    'dipoles_nuc',
    'btensor',
    'ctensor',
    'magnetizability',
    'dtensor',
    'dstensor',
    'tdspec_wavelength',
    'tdspec_wavenumber',
    'tdspec_opa',
    'tdspec_rrs',
    'tdspec_tpa',
    'tdspec_rhrs',
    'tdspec_drsfg',
    'atomic_charges',
    'atomic_dipoles',
    'atomic_quadrupoles',
    'density',
    'fde',
    'hirshfeld_induced_dipoles_loc',
    'hirshfeld_induced_dipoles_nonloc',
    'hirshfeld_induced_dipoles_tot',
    'hirshfeld_induced_charges',
    'hirsh_pol',
    'beta_ddq',
    'beta_dqd',
    'beta_qdd',
    'beta_dqq',
    'beta_qdq',
    'beta_qqd',
    'beta_qqq',
    'project', # added to save time for TERS - Pengchong, Nov. 2016
    'basis',
)

attrset = frozenset(attrlist)

privates = frozenset([
    '_i',
    '_ct',
    '_diagonalized',
    '_excite_index',
    '_raman',
    '_input_order',
    '_privates',
    '_attrset',
    '_attrlist',
    '_managed_methods',
    '_computer_names',
    '_abort',
])

managed_methods = (
    'allcoords',
    'allatoms',
    'allelements',
    'nallatoms',
    'nallelements',
    'masses',
    'dim_masses',
    'allmasses',
    'molecular_mass',
    'dim_mass',
    'system_mass',
    'bonds',
    'dim_bonds',
    'polarizability'
)

computer_names = (
     'Helium',
     'Lithium',
     'Beryllium',
     'Boron',
     'Carbon',
     'Nitrogen',
     'Oxygen',
     'Fluorine',
     'CyberStar',
     'LionXJ',
     'LionXF',
     'LionXC',
     'LionXK',
     'LionXI',
     'Hammer',
)
