from __future__ import print_function, division
from numpy import array, append

def collect_polarizability(self, f, indices):
    '''Drive the collection of polarizabilities from different methods.'''

    if 'DFT POL' in indices:
        __dft_polarizability(self, f, indices)
    elif 'CCSD POL' in indices:
        __ccsd_polarizability(self, f, indices)
    # Add other methods here
    #elif 'OTHER POL' in indices:
    else:
        self._raise_or_pass('No polarizability block was found')

def __dft_polarizability(self, f, indices):
    '''Collect the DFT polarizabilities.'''

    # For frequency-dependent polarizability calculations, damping
    # (lifetime) in the calculation causes the polarizability to
    # be complex.
    ar = indices['DFT POL']
    self.e_frequencies = array([], dtype=float)
    self.qm_pol = array([], dtype=float)
    # Loop over all possible locations for this header
    for ix in ar:
        ln = f[ix+1].split()
        self.e_frequencies = append(self.e_frequencies, float(ln[2]))
        if 'DAMPING' in self.subkey:
            # Real
            s = ix + 5
            e = ix + 8 
            r = array([x.split()[1:4] for x in f[s:e]], dtype=float)
            # Imaginary
            s = ix + 17
            e = ix + 20 
            # Need a try statement here because it is possible for the
            # numbers to overlap.
            try:
                i = array([x.split()[1:4] for x in f[s:e]], dtype=float)
            except ValueError:
                i = array([[x[5:15],x[15:24],x[24:33]] for x in f[s:e]], dtype=float)
            pol = r + i*1j
        else:
            pol = f[ix+5:ix+8]
            pol = array([x.split()[1:4] for x in pol], dtype=float)
        try: 
            self.qm_pol = append(self.qm_pol, array([pol]), axis=0)
        except ValueError:
            self.qm_pol = array([pol])
    self.npol = len(self.e_frequencies)

def __ccsd_polarizability(self, f, indices):
    '''Collect the CCSD polarizabilities.'''

    ar = indices['CCSD POL']
    self.e_frequencies = array([], dtype=float)
    self.qm_pol = array([], dtype=float)
    # Loop over all possible locations for this header (unlikely that
    # more than one is present).
    for ix in ar:
        ln = f[ix+1].split() 
        self.e_frequencies = append(self.e_frequencies, float(ln[2]))
        pol = f[ix+7:ix+10]
        pol = array([x.split()[1:4] for x in pol], dtype=float)
        try:
            self.qm_pol = append(self.qm_pol, array([pol]), axis=0)
        except ValueError:
            self.qm_pol = array([pol])
    self.npol = len(self.e_frequencies)
