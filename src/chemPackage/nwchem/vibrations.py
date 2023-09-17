from __future__ import print_function, division
from numpy import array, where, append, argsort, row_stack
from ..constants import ANGSTROM2BOHR as A2B
from ..constants import BOHR2ANGSTROM as B2A
import numpy

def collect_frequencies(self, f, indices):
    '''Collect frequencies and IR intensities.'''

    # Function to find the block separators
    def vib_freq(line):
        ln = line.split()
        if not ln:
            return False
        elif 'P.Frequency' in ln:
            return True
        else:
            return False

    # Find the limits of the frequencies block
    if 'NMODES START' in indices:
        s = indices['NMODES START']
        e = indices['NMODES END']

        # Locate where each header is
        ar = [i for i, x in enumerate(f[s:e], s) if vib_freq(x)]

        # Collect the normal modes and vibrational frequencies
        self.v_frequencies = array([], dtype=float)
        self.normal_modes = array([], dtype=float)
        for s in ar:

            # Collect the vibrational frequencies for this set
            ln = array(f[s].split()[1:], dtype=float)
            self.v_frequencies = append(self.v_frequencies, ln)

            # Collect the mode for this set
            s += 2
            e = s + self.natoms * 3
            # Place into array holding only the modes (remove first column)
            modes = array([x.split() for x in f[s:e]], dtype=float)[:,1:].T
            # Reshape so that modes are 2D, not 1D
            modes = modes.reshape((len(ln), self.natoms, 3))
            # Add to the normal modes array
            try:
                self.normal_modes = append(self.normal_modes, modes, axis=0)
            except ValueError:
                self.normal_modes = modes

        # Finishing touches
        self.nmodes = len(self.v_frequencies)
        # Simultaneously normalize in Bohr then convert to Angstroms
        for i in range(self.nmodes):
            self.normal_modes[i] *= B2A / numpy.linalg.norm(self.normal_modes[i].flatten())
    else:
        self._raise_or_pass('Error locating normal modes')

    # Collect the IR intensities
    if 'IR' in indices:
        s = indices['IR']
        e = next(i for i, x in enumerate(f[s+1:], s+1) 
                                                    if '----------------' in x)
        self.IR = array([x.split()[5] for x in f[s:e]], dtype=float)
    else:
        self._raise_or_pass('Error locating IR intensities')
