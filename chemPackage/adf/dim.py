from __future__ import print_function, division
from copy import deepcopy
from numpy import array, append, asarray
from numpy import fastCopyAndTranspose as fcat

def collect_dim(self, f, indices):
    '''Collect all DIM properties.'''

    ########################
    # COLLECT DIM CORDINATES
    ########################

    # The DIM coordinates are laid out similarly to the QM coordinates,
    # so we are able to use the same regular expression to collect
    # them.
    if 'DIM COORDINATES' in indices:
        s = indices['DIM COORDINATES']
        e = next(i for i, x in enumerate(f[s:], s) if not x.strip())
        # Elements 2, 3 and 4 are the X, Y, and Z coordinates
        self.dim_coordinates = array([x.split()[2:5] for x in f[s:e]],
                                                           dtype=float)
        self.dim_atoms = array([x.split()[1] for x in f[s:e]])
        self.dim_elements = set(self.dim_atoms)
        self.ndim_elements = len(self.dim_elements)
        self.ndim = len(self.dim_atoms)
        from ..constants import BOHR2ANGSTROM
        self.dim_coordinates = BOHR2ANGSTROM(self.dim_coordinates)
    else:
        self._raise_or_pass('Error locating DIM coordinates')

    #########################################
    # COLLECT DIM INDUCED CHARGES AND DIPOLES
    #########################################

    # The DIM induced charges and dipoles are listed in a table with
    # the dipole vector first and then the charge.  
    if 'DIM DIPOLES AND CHARGES' in indices:
        start = indices['DIM DIPOLES AND CHARGES']
    elif 'DIM DIPOLES' in indices:
        start = indices['DIM DIPOLES']
    elif 'OLD DIM D&C' in indices:
        start = indices['OLD DIM D&C']
    else:
        start = None
    # Read in the dipoles and charges if they are in the output file
    if start:
        self.dim_dipoles = {}
        if 'CPIM' in self.subkey: self.dim_charges = {}

        # Locate all locations of the induced charges and dipoles
        for st in start:

            # Start and stop limits
            s = next(i for i, x in enumerate(f[st:], st) if 'ATOM' in x) + 3
            e = next(i for i, x in enumerate(f[st:], st) if '========' in x)

            # Complex boolean
            cmplx = ('FD' in self.calctype and len(f[st+6].split()) == 6)
            #cmplx = ('FD' in self.calctype and len(f[st+7].split()) == 5)
            # Complex.
            if cmplx:
                # Elements 2, 3, and 4 on every other line are the real dipoles
                r = array([x.split()[2:5] for x in f[s:e-1:2]], dtype=float)
                # Elements 0, 1, and 2 on every other line are the imag dipoles
                i = array([x.split()[0:3] for x in f[s+1:e:2]], dtype=float)
                dpls = r + i*1j
                # Charge not printed for PIM and DRF.
                if 'CPIM' in self.subkey:
                    # Element 5 on every other line is the real charge.
                    r = array([x.split()[5] for x in f[s:e-1:2]], dtype=float)
                    # Element 3 on every other line is the imag charge.
                    i = array([x.split()[3] for x in f[s+1:e:2]], dtype=float)
                    chrgs = r + i*1j

            # Real
            else:
                #print(array([x.split()[2:5] for x in f[s:e]], dtype=float))
                dpls = array([x.split()[2:5] for x in f[s:e]], dtype=float)
                if 'CPIM' in self.subkey:
                    chrgs = array([x.split()[5] for x in f[s:e]], dtype=float)

            # Assign to proper location based on description after header
            if '---------------' in f[st+1]:
                self.dim_dipoles['ground image'] = dpls
                if 'CPIM' in self.subkey:
                    self.dim_charges['ground image'] = chrgs
            elif not cmplx:
                head = 'static scattered'
                # Store the different direction sequentially for each frequency
                dir = { 'X': 0, 'Y': 1, 'Z': 2 }[f[st+1].split()[0]]
                try:
                    self.dim_dipoles[head][dir] = dpls
                    if 'CPIM' in self.subkey: self.dim_charges[head][dir] = chrgs
                except KeyError:
                    self.dim_dipoles[head] = [None, None, None]
                    self.dim_dipoles[head][dir] = dpls
                    if 'CPIM' in self.subkey:
                        self.dim_charges[head] = [None, None, None]
                        self.dim_charges[head][dir] = chrgs
            else:
                head = 'FD scattered' if 'Local' in f[st+1] else 'ES image'
                # Store the different direction sequentially for each frequency

                dir = { 'X': 0, 'Y': 1, 'Z': 2 }[f[st+1].split()[0]]
                try:
                    self.dim_dipoles[head][dir].append(dpls)
                    if 'CPIM' in self.subkey:
                        self.dim_charges[head][dir].append(chrgs)
                except KeyError:
                    self.dim_dipoles[head] = [[], [], []]
                    self.dim_dipoles[head][dir].append(dpls)
                    if 'CPIM' in self.subkey:
                        self.dim_charges[head] = [[], [], []]
                        self.dim_charges[head][dir].append(chrgs)

        # Reorder so that it goes head:freq:dir:atom(:xyz) instead
        # of head:dir:freq:atom(:xyz)
        for head in self.dim_dipoles:
            if head == 'ground image':
                continue
            elif head == 'static scattered':
                self.dim_dipoles[head] = array(self.dim_dipoles[head])
                if 'CPIM' in self.subkey:
                    self.dim_charges[head] = array(self.dim_charges[head])
            else:
                try:
                    self.dim_dipoles[head] = array(self.dim_dipoles[head]).swapaxes(0,1)
                except ValueError:
                # Temporarily disabled ValueError raise; 
                # if calculated only one direction in AOResponse, other dir have empty list, 
                # which leads to axes swap ValueErrors. We can ignore this case. --Pengchong, Jan 2018
                    continue;
                if 'CPIM' in self.subkey:
                    self.dim_charges[head] = array(
                       self.dim_charges[head]).swapaxes(0,1)

    ####################################
    # COLLECT DIM/QM TOTAL DIPOLE MOMENT
    ####################################

    # If this is a DIMQM run, then we want to collect the DIM total
    # dipole moment.  It is after the results header just after the line
    # discussing the total dipole moment.
    if 'DIM DIPOLE MOMENT' in indices:
        s = indices['DIM DIPOLE MOMENT']
        sl = ' Total induced dipole moment in DIM system :'
        ix =  next(i for i, x in enumerate(f[s:], s) if x == sl)
        self.dim_dipole_tot = array(f[ix+3].split()[1:4], dtype=float)
    else:
        self._raise_or_pass('Error locating DIM total dipole moment')

    ###################################
    # COLLECT DIM/QM INTERACTION ENERGY
    ###################################

    # If this is a DIMQM run, then we want to collect the DIM/QM
    # interaction energy.  It on the line with 'Total' after the header
    # given below
    if 'DIM/QM ENERGY' in indices:
        s = indices['DIM/QM ENERGY']
        ix = next(i for i, x in enumerate(f[s:], s) if '  Total' in x)
        self.dimqm_energy = float(f[ix].split()[2])
    else:
        self._raise_or_pass('Error locating DIM/QM energy')
        
    ###########################
    # COLLECT DIM SYSTEM ENERGY
    ###########################

    # If this is a DIMQM run, then we want to collect the DIM system
    # energy.  If it on the line with 'Total' after the header given below
    if 'DIM ENERGY' in indices:
        s = indices['DIM ENERGY']
        ix = next(i for i, x in enumerate(f[s:], s) if '  Total' in x)
        self.dim_energy = float(f[ix].split()[2])
        
    ####################################
    # COLLECT DIM POLARIZABILITY TENSORS
    ####################################

    # The polarizability tensor for the DIM system is collected.
    # It can be real or complex, depending on the calculation.
    if 'POLARIZABILITY' in self.calctype:
        if 'FD' in self.calctype:
            if 'DIM FD REAL POLARIZABILITY' in indices:
                arr = indices['DIM FD REAL POLARIZABILITY']
            else:
                self._raise_or_pass('Error locating DIM real polarizability')
                arr = []
            if 'DIM FD IMAG POLARIZABILITY' in indices:
                # Make Major axis the frequencies
                arr = fcat(array([arr,indices['DIM FD IMAG POLARIZABILITY']]))
            else:
                self._raise_or_pass('Error locating DIM imag polarizability')
                arr = []
        else:
            if 'DIM STATIC POLARIZABILITY' in indices:
                arr = indices['DIM STATIC POLARIZABILITY']
            else:
                self._raise_or_pass('Error locating DIM static polarizability')
                arr = []

        for ar in arr:
            if 'FD' in self.calctype:
                # Collect complex isolated DIM tensor
                s = ar[0] + 3
                e = ar[0] + 6
                r = array([[x.split()[1:4] for x in f[s:e]]], dtype=float)
                s = ar[1] + 3
                e = ar[1] + 6
                i = array([[x.split()[1:4] for x in f[s:e]]], dtype=float)
                try:
                    self.dim_pol = append(self.dim_pol, r+i*1j, axis=0)
                except ValueError:
                    self.dim_pol = r+i*1j
            else:
                # Collect real isolated DIM tensor
                s = ar + 3
                e = ar + 6
                r = array([[x.split()[1:4] for x in f[s:e]]], dtype=float)
                try:
                    self.dim_pol = append(self.dim_pol, r, axis=0)
                except ValueError:
                    self.dim_pol = r
