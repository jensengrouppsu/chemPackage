from __future__ import print_function, division
from numpy import array, append, reshape, zeros

def collect_optical_rotation(self, f, indices):
    '''Drive the collection of optical rotation tensors from different methods.'''

    if 'DFT POL' in indices:
        __dft_optical_rotation(self, f, indices)
    # Add other methods here
    #elif 'OTHER POL' in indices:
    else:
        self._raise_or_pass('No optical rotation block was found')

def __dft_optical_rotation(self, f, indices):
    '''Collect the DFT optical rotation tensors.'''

    # The incident frequency is already collected during the
    # polarizability calculation of the response.  Also, the
    # number of tensors must match the number of polarizability
    # tensors, so we have that already also.
    ar = indices['DFT OPTICAL ROTATION']
    self.ord = array([], dtype=float)
    # Loop over all possible locations for this header
    for ix in ar:
        if 'DAMPING' in self.subkey:
            # Real
            s = ix + 4
            e = ix + 7 
            try:
                r = array([x.split()[1:4] for x in f[s:e]], dtype=float)
            except ValueError: #No space between columns due to overflow
                import re
                re1='.*?'
                re2='([-]?\\d*\\.\\d+)'
                re3='\\s*'
                rg = re.compile(re1+re2+re3+re2+re3+re2,re.IGNORECASE|re.DOTALL)
                r = array([], dtype=float)
                for x in f[s:e]:
                    m = rg.search(x)
                    try:
                        t = array([m.group(1), m.group(2), m.group(3)], dtype=float)
                        r = append(r, t)
                    except AttributeError:
                        #Overflow of values in Fortran (printed ******)
                        r = zeros((3,3))
                        break
                r = reshape(r, (3,3))
            # Imaginary
            s = ix + 11
            e = ix + 14 
            try:
                i = array([x.split()[1:4] for x in f[s:e]], dtype=float)
            except ValueError:
                import re
                re1='.*?'
                re2='([-]?\\d*\\.\\d+)'
                re3='\\s*'
                rg = re.compile(re1+re2+re3+re2+re3+re2,re.IGNORECASE|re.DOTALL)
                i = array([], dtype=float)
                for x in f[s:e]:
                    m = rg.search(x)
                    t = array([m.group(1), m.group(2), m.group(3)], dtype=float)
                    i = append(i, t)
                i = reshape(i, (3,3))
            ord = r + i*1j
        else:
            ord = f[ix+4:ix+7]
            try:
                ord = array([x.split()[1:4] for x in ord], dtype=float)
            except ValueError:
                ord = array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], dtype=float)
        try: 
            self.ord = append(self.ord, array([ord]), axis=0)
        except ValueError:
            self.ord = array([ord])
