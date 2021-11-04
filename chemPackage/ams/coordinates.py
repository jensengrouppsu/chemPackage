from numpy import where, array, append
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
    



