from __future__ import print_function, division
from numpy import where, array, append
import re

def collect_geometry(self, f, indices):
    '''Collect the initial geometry of the molecule for the calculation.'''

    # Find number of atoms
    if 'ATOMS NUMBER' in indices:
        ix = indices['ATOMS NUMBER']
        self.natoms = int(f[ix].strip())
    else:
        self.__raise_or_pass('Could not find number of atoms')
    self.natoms = int(f[ix].strip())

    # Locate coordinates block
    if 'GEOMETRY' in indices:
        s = indices['GEOMETRY'][0]
        e = s + self.natoms
        # Collect the coordinate data
        self.coordinates = array([x.split()[3:6] for x in f[s:e]], dtype=float)
        self.atoms = array([x.split()[1] for x in f[s:e]])
        # Set the elements in the system
        self.elements = set(self.atoms)
        self.nelements = len(self.elements)
    else:
        self._raise_or_pass('Could not find initial coordinates.')


def collect_optimized_geometry(self, f, indices):
    '''Collect geometry optimization coordinates, both final and initial.'''

    # Copy old coordinates to initial geometry
    self.initial_geometry = self.coordinates.copy()

    # Find new coordinates
    if 'GEOMETRY' in indices:
        s = indices['GEOMETRY'][-1]
        e = s + self.natoms
        self.coordinates = array([x.split()[3:6] for x in f[s:e]], dtype=float)
    else:
        self._raise_or_pass('Could not find optimized coordinates.')


def geometry_block(self, s, e, f):
    '''Collect the coordinates and atoms from the GEOMETRY block.
    The coordinates are stored as floats. This will be overwritten unless
    this is an input file.
    '''

    # Dictionary of compiled regular expressions for faster excecution
    coordinate = re.compile(
                            r'''
                             \s*           # Leading whitespace
                             ([a-z]?[A-Z]?[a-z]?[0-9]?) # The element or charge designation
                             \s+           # Whitespace
                             ([-]?[0-9.]+([eE][-+]?[0-9]+)?)    # X coord (posibly in scientific notation)
                             \s+           # etc...
                             ([-]?[0-9.]+([eE][-+]?[0-9]+)?)
                             \s+
                             ([-]?[0-9.]+([eE][-+]?[0-9]+)?)
                            ''', re.VERBOSE)
    #[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?
    # Put the coordinates and atoms from ATOMS block into the correct place
    self.coordinates = []
    self.atoms = []
    for ln in f[s:e]:
        if 'SYMMETRY' in ln.upper() or not ln.strip(): continue
        m = coordinate.match(ln)
        try:
            self.atoms = append(self.atoms, [m.group(1)])
            self.coordinates = append(self.coordinates,
                                          [m.group(2), m.group(4), m.group(6)])
        except (IndexError, AttributeError):
            self._raise_or_pass('Error found in GEOMETRY block')
    self.natoms = len(self.atoms)
    self.coordinates = array(self.coordinates, dtype=float).reshape(-1,3)
    # Set the elements in the system
    self.elements = set(self.atoms)
    self.nelements = len(self.elements)

def collect_ground_state_gradient(self, f, indices):
    '''Collect the ground state gradient for a DFT calculation.'''

    # Initialize the gradient dictionary and location dictionary
    self.gs_gradient = {}
    gsg_locations = {}

    # Table of things to search for
    table = ('nuclear repulsion gradient',
             'weighted density gradient',
             'kinetic energy gradient',
             '2-electron gradient',
             'DFT CD+XC gradient',)

    # Find the gradients
    if 'DFT GRAD START' in indices:
        s = indices['DFT GRAD START'][-1]
        if len(indices['DFT GRAD START']) != len(indices['DFT GRAD END']):
            e = indices['DFT GRAD END'][-2] 
        else:
            e = indices['DFT GRAD END'][-1] 
        count = s
        for ln in f[s:e]:
            # Find locations of contributions to the gradient (if they're printed)
            for item in table:
                if item in ln:
                    gsg_locations[item] = count
            # Collect the total ground state gradient
            if 'DFT ENERGY GRADIENTS' in ln:
                self.gs_gradient['total DFT gradient'] = (
                    array([x.split()[5:8] for x in f[count+4:count+4+self.natoms]], dtype=float))
            count += 1
        # Collect the contributions to the gradient
        for item in gsg_locations.keys():
            s = gsg_locations[item] + 1
            e = gsg_locations[item] + 1 + self.natoms
            self.gs_gradient[item] = array([x.split()[0:3] for x in f[s:e]], dtype=float)
    else:
        self._raise_or_pass('Could not find ground state gradient.')
