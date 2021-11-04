from numpy import array, where, append, arange, zeros
from numpy import row_stack, column_stack, argsort
from chemPackage import collect

def collect_polarizability(self, f, indices):
    '''Drive the collection of polarizabilities from different methods.'''
    __aoresponse(self, f, indices)



def __aoresponse(self, f, indices):
    '''Collect AORESPONSE polarizability.'''

    # This section collections the AORESPONSE routine tensors. Each
    # polarizability is identified by the frequency in Hartrees that
    # it was calculated at.  
    #
    # AORESPONSE can report imaginary polarizability tensors.  If the
    # imaginary tensor is reported, then it is added to the real to
    # make it complex.
    if 'AORESPONSE' in indices:
        ar = indices['AORESPONSE']
        self.e_frequencies = array([], dtype=float)
        for ix in ar:
            # Collect frequency
            ln = f[ix+1].split()
            self.e_frequencies = append(self.e_frequencies, float(ln[2]))
            # Collect tensor.  Different for real and complex.
            if f[ix+14] == "        IMAGINARY POLARIZABILITY":
           #if 'LIFETIME' in self.subkey:
                s = ix + 8
                e = ix + 11
                r = array([[x.split() for x in f[s:e]]], dtype=float)
                try:
                    s = ix + 17
                    e = ix + 20
                    i = array([[x.split() for x in f[s:e]]], dtype=float)
                except ValueError:
                    s = ix + 16
                    e = ix + 19
                    i = array([[x.split() for x in f[s:e]]], dtype=float)
                pol = r + i*1j
            else:
                s = ix + 8
                e = ix + 11
                pol = array([[x.split() for x in f[s:e]]], dtype=float)
            # Add this tensor to the array.
            try: 
                self.qm_pol= append(self.qm_pol, pol, axis=0)
            except ValueError:
                self.qm_pol = pol
        self.npol = len(self.e_frequencies)
