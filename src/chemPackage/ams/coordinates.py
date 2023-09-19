import numpy as np
import re

def collect_geometry(self, f, indices):
    '''Collect the geometry of the molecule used for the calculation.
    This is not used for a geometry optimization.
    '''

    # The coordinates are printed at the top of the output file
    if 'INITIAL GEOMETRY' in indices or 'INITIAL GEOMETRY COSMO' in indices:    
        # For some reason, the initial geometry in COSMO is printed differently
        if 'INITIAL GEOMETRY' in indices:
            s = indices['INITIAL GEOMETRY']
        else:
            s = indices['INITIAL GEOMETRY COSMO']

        # find number of atoms in system
        natoms = 0
        while len(f[s+natoms].split()) == 8:
            natoms = natoms+ 1
        self.natoms = natoms
        # Go to end of block
        try:
            e = s + self.natoms
        except TypeError:
            self._raise_or_pass('Error collecting initial geometry block')
            return
        # Elements 2, 3, and 4 are the X, Y, and Z coordinates
        self.coordinates = array([x.split()[2:5] for x in f[s:e]],dtype=float)
        self.atoms = array([x.split()[1] for x in f[s:e]])
        # Store these as initial geometry as well
        self.initial_geometry = self.coordinates.copy()
    else:
        self._raise_or_pass('Error locating initial geometry block')

def collect_optimized_geometry(self, f, indices):
    '''Collect geometry optimization coordinates, both final and initial.'''

    # This section will collect the molecular coordinates from the a
    # geometry optimization.  Because the geometry optimization tables print
    # out more digits than the initial geometry table, we will collect the
    # coordinates from the first and last optimization table, which gives
    # the initial and final coordinates, respectively.  
    if 'OPTIMIZED GEOMETRY' in indices:
        # Initial geometry
        s = indices['OPTIMIZED GEOMETRY'][0]
        try:
            e = s + self.natoms
        except TypeError:
            self._raise_or_pass('Error collecting optimized geometry block')
            return
        # Elements 5, 6 and 7 are the X, Y, and Z coordinates
        self.initial_geometry = array([x.split()[2:5] for x in f[s:e]], dtype=float)

        # Final geometry
        s = indices['OPTIMIZED GEOMETRY'][-1]
        e = s + self.natoms
        # Elements 5, 6 and 7 are the X, Y, and Z coordinates
        self.coordinates = array([x.split()[2:5] for x in f[s:e]],dtype=float)
    else:
        self._raise_or_pass('Error locating optimized geometry blocks')
    
def atom_block(self, s, e, search, f, indices):
    '''Collect the coordinates and atoms from the ATOMS block.
    The coordinates are stored as floats, after the possible substitution
    from a GEOVAR block.  This will be overwritten unless this is an input
    file.
    '''

    # Dictionary of compiled regular expressions for faster excecution
    coordinate = re.compile(
                            r'''
                             \s*           # Leading whitespace
                             \d*           # Atomic ## (optional)
                             .*            # Separator
                             ([A-Z][a-z]?(.Gh)?) # The element -- allow for ghost atoms
                             \s+           # Whitespace
                             ([-0-9.]+([eE][-+]?[0-9]+)?)   # X coord or variable
                             \s+           # etc...
                             ([-0-9.]+([eE][-+]?[0-9]+)?)
                             \s+
                             ([-0-9.]+([eE][-+]?[0-9]+)?)
                             (\s*f=.*)?    # maybe geometry is from fragment files
                            ''', re.VERBOSE)

    # Put the coordinates and atoms from ATOMS block into the correct place
    self.coordinates = []
    self.atoms = []
    for ln in f[s:e]: #xing
        m = coordinate.match(ln)
        try:
            self.atoms = np.append(self.atoms, [m.group(1)])
            self.coordinates = np.append(self.coordinates,
                                          [m.group(3), m.group(5), m.group(7)])
        except (IndexError, AttributeError):
            self._raise_or_pass('Error found in ATOMS block: ' + ln)
    self.natoms = len(self.atoms)
    self.coordinates = self.coordinates.reshape(-1,3)
    # Set the elements in the system
    self.elements = set(self.atoms)
    self.nelements = len(self.elements)

    # Make geovar subsitutions to the coordinates now.
    s = indices['INPUT START']
    e = indices['INPUT END']
    sx = next((i for i,x in enumerate(search[s:e],s) if x == 'GEOVAR'), -1) + 1
    if sx:
        ex = next(i for i, x in enumerate(search[s:e], s) if x == 'END')
        for g in f[sx:ex]:
            ln = g.split()
            try:
                ix = np.where(self.coordinates == ln[0])
                self.coordinates[ix[0][0]][ix[1][0]] = ln[1]
            except IndexError:
                self._raise_or_pass('Error collecting GEOVAR block: ' + g)

    # Make coordinates into floats
    self.coordinates = np.array(self.coordinates, dtype=float)

