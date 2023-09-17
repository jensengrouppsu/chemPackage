from __future__ import print_function, division
from numpy import array, append

def collect_occ_vir_mos(self, f, indices):
    '''Collect the number of occupied and virtual MOs.'''

    try:
        __collect_occ_vir(self, f, indices)
    except IndexError:
        self._raise_or_pass('Error collecting number of occ./vir. MOs')


def __collect_occ_vir(self, f, indices):
    '''Collect numbers of MOs from the general information block.'''

    # Collection is identical regardless of whether a restricted or 
    # unrestricted calculation is done. 

    start = [indices['MO NUMBER']]

    # Allocate the occupied and virtual orbital storing arrays.
    self.nocc = array([], dtype=int)
    self.nvir = array([], dtype=int)

    noccbeta = 0

    # Loop for collection
    for n, s in enumerate(start):
        
        # Find total number of MOs and occupied MOs
        for i, x in enumerate(f[s:], s):
            # Number of alpha MOs
            if 'Alpha electrons' in x:
                # For DFT
                noccalpha = int(x.split()[3])
            elif 'closed shells' in x:
                # For RHF and ROHF
                noccalpha = int(x.split()[3])
            elif 'alpha electrons' in x:
                # For UHF 
                noccalpha = int(x.split()[3])
            # Number of beta MOs
            if 'Beta electrons' in x:
                # For DFT
                noccbeta = int(x.split()[3])
            elif 'beta electrons' in x:
                # For UHF
                noccbeta = int(x.split()[3])  
            # Number of MOs
            if 'AO basis' in x:
                # For DFT
                nmo = int(x.split()[6])
            elif 'functions' in x:
                # For RHF/ROHF/UHF
                nmo = int(x.split()[2])
            if 'XC Information' in x:
                break
            elif 'Summary' in x:
                # For RHF, ROHF, UHF 
                break
        if noccbeta == 0:
            # For RHF
            noccbeta = noccalpha

        # Find the number of virtual MOs
        nviralpha = nmo - noccalpha
        nvirbeta = nmo - noccbeta
 
    # Store MO numbers.  For restricted calculations, the numbers of
    # occupied and virtual orbitals are identical for alpha/beta spin.
    # For unrestricted calculations, the numbers may be the same or
    # different. 
    self.nocc = append(self.nocc, noccalpha)
    self.nocc = append(self.nocc, noccbeta)
    self.nvir = append(self.nvir, nviralpha)
    self.nvir = append(self.nvir, nvirbeta)

    # Store the numbers of alpha/beta occupied and virtual orbitals as
    # a tuple.
    self.nocc = tuple(self.nocc)
    self.nvir = tuple(self.nvir)

