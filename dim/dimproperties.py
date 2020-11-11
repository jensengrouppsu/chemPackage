from __future__ import print_function, division
from numpy import array, where, zeros, append, argsort
from numpy import fastCopyAndTranspose as fcat
from ..constants import BOHR2ANGSTROM

def collect_dim(self, f, indices):
    '''Collect all DIM properties.'''

    ########################
    # COLLECT DIM CORDINATES
    ########################
    if 'COORDS' in indices:
        s = indices['COORDS']
        e = next(i for i, x in enumerate(f[s:], s) if not x.strip())
        # Elements 2, 3 and 4 are the X, Y, and Z coordinates
        self.dim_coordinates = array([x.split()[2:5] for x in f[s:e]],
                                                           dtype=float)
        self.dim_atoms = array([x.split()[1] for x in f[s:e]])
        self.dim_elements = set(self.dim_atoms)
        self.ndim_elements = len(self.dim_elements)
        self.ndim = len(self.dim_atoms)
        self.dim_coordinates = BOHR2ANGSTROM(self.dim_coordinates)
        # Intialize natoms (QM atoms) to 0
        self.natoms = 0
        self.nelements = 0
        # Collect coordination radii from element 5
        if 'COORDDEPEND' in self.key:
            self.dim_cnradii = array([x.split()[5] for x in f[s:e]], dtype=float)
            self.dim_cnradii = BOHR2ANGSTROM(self.dim_cnradii)
        if "COARSEGRAIN" in self.key:
            self.dim_cgradii = array([x.split()[5] for x in f[s:e]], dtype=float)
            self.dim_cgradii = BOHR2ANGSTROM(self.dim_cgradii)
    # Try to collect from the xyz file
    '''else:
        try:
            from ..__init__ import collect
            tmp = collect(self.key['XYZ'][0])
            tmp.make_dim_atoms()
            # Shorter to directly copy over xyz info than list off what not to copy
            from copy import deepcopy
            for k in ['ndim', 'ndim_elements', 'dim_atoms', 'dim_elements', 'dim_coordinates']:
                setattr(self, k, deepcopy(getattr(tmp, k)))
        except IOError:
            if 'RUNTIME' in indices:
                ix  = indices['RUNTIME'] + 13
                self.ndim = int(f[ix].split()[3])
                # Intialize natoms (QM atoms) to 0
                self.natoms = 0
                self.nelements = 0'''

    ###################################
    # DETERMINE CALCULATION FREQUENCIES
    ###################################

    # Define the search string
    if 'FREQUENCY' in indices:
        ar = indices['FREQUENCY']
        self.e_frequencies = array([], dtype=float)
        for ix in ar:
            # Collect the frequency
            ln = f[ix].split()
            self.e_frequencies = append(self.e_frequencies, float(ln[1]))
        # Store the number of frequencies
        self.npol = len(self.e_frequencies)
    else:
        self._raise_or_pass('Error locating DIM frequencies')

    #########################################
    # COLLECT DIM INDUCED CHARGES AND DIPOLES
    #########################################

    if 'CPIMDIP' in indices:
        start = indices['CPIMDIP']
    elif 'PIMDIP' in indices:
        start = indices['PIMDIP']
    else:
        start = []

    if start:
        self.dim_dipoles = {}
        if 'CPIM' in self.subkey: self.dim_charges = {}

        # Locate all locations of the induced charges and dipoles
        for st in start:

            # Start and stop limits
            s = st + 3
            e = next(i for i, x in enumerate(f[s:], s) if not x.strip())

            # Complex.
            if 'FD' in self.calctype:
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

                # Assign to proper location based on description after header
                head = 'FD scattered'
                # Store the different direction sequentially for each frequency
                dir = { 'X': 0, 'Y': 1, 'Z': 2 }[f[st-3].split()[-3]]
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

            # Real
            else:
                dpls = array([x.split()[2:5] for x in f[s:e]], dtype=float)
                if 'CPIM' in self.subkey:
                    chrgs = array([x.split()[5] for x in f[s:e]], dtype=float)

                # Assign to proper location based on description after header
                head = 'static scattered'
                # Store the different direction sequentially for each frequency
                dir = { 'X': 0, 'Y': 1, 'Z': 2 }[f[st-3].split()[-3]]

                try:
                    self.dim_dipoles[head][dir] = dpls
                    if 'CPIM' in self.subkey: self.dim_charges[head][dir] = chrgs
                except KeyError:
                    self.dim_dipoles[head] = [None, None, None]
                    self.dim_dipoles[head][dir] = dpls
                    if 'CPIM' in self.subkey:
                        self.dim_charges[head] = [None, None, None]
                        self.dim_charges[head][dir] = chrgs

        # Reorder so that it goes head:freq:dir:atom(:xyz) instead
        # of head:dir:freq:atom(:xyz)
        for head in self.dim_dipoles:
            if head == 'static scattered':
                self.dim_dipoles[head] = array(self.dim_dipoles[head])
                if 'CPIM' in self.subkey:
                    self.dim_charges[head] = array(self.dim_charges[head])
            else:
                self.dim_dipoles[head] = array(
                    self.dim_dipoles[head]).swapaxes(0,1)
                if 'CPIM' in self.subkey:
                    self.dim_charges[head] = array(
                       self.dim_charges[head]).swapaxes(0,1)

    # If this table was not found but we were supposed to, alert.
    else:
        if 'ATMDIP' in self.subkey:
            self._raise_or_pass('Could not find DIM atomic dipoles.')
        else:
            pass


    ####################################
    # COLLECT DIM POLARIZABILITY TENSORS
    ####################################

    # The polarizability tensor for the DIM system is collected.
    # It can be real or complex, depending on the calculation.
    if 'POL' in indices:
        ar = indices['POL']
        for ix in ar:
            if 'FD' in self.calctype:
                # Collect complex isolated DIM tensor
                s = ix + 5
                e = ix + 8
                r = array([[x.split()[1:4] for x in f[s:e]]], dtype=float)
                s = ix + 11
                e = ix + 14
                i = array([[x.split()[1:4] for x in f[s:e]]], dtype=float)
                try:
                    self.dim_pol = append(self.dim_pol, r+i*1j, axis=0)
                except ValueError:
                    self.dim_pol = r+i*1j
            else:
                # Collect real isolated DIM tensor
                s = ix + 4
                e = ix + 7
                r = array([[x.split()[1:4] for x in f[s:e]]], dtype=float)
                try:
                    self.dim_pol = append(self.dim_pol, r, axis=0)
                except ValueError:
                    self.dim_pol = r

    # If this table was not found but we were supposed to, alert.
    else:
        if 'POLARIZABILITY' in self.calctype:
            self._raise_or_pass('Could not find DIM polarizabilities')
        else:
            pass

    #########################################
    # COLLECT EFFICIENCIES AND CROSS SECTIONS
    #########################################

    if 'EFF' in indices:
        ar = indices['EFF']
        for ix in ar:
            ln = f[ix].split()
            ln = [ float(ln[0]), float(ln[1]), float(ln[2]),
                   float(ln[3]), float(ln[4]), float(ln[5]) ]
            # Need to store efficiencies as an array of lists, since the
            # reordering unintentionally removes all array elements beyond
            # the element indexed by the number of frequencies.
            try:
                self.dim_efficiencies   = append(self.dim_efficiencies,
                                                 [ ln[0:3] ], axis=0)
                self.dim_cross_sections = append(self.dim_cross_sections,
                                                 [ ln[3:6] ], axis=0)
            except ValueError:
                self.dim_efficiencies   = array([ ln[0:3] ], dtype=float)
                self.dim_cross_sections = array([ ln[3:6] ], dtype=float)

    # Put the frequencies in ascending order
    indx = argsort(self.e_frequencies)
    # Sort the other properties accordingly
    if self.dim_dipoles is not None:
        for head in self.dim_dipoles:
            if head is not 'static scattered':
                self.dim_dipoles[head] = self.dim_dipoles[head][indx]
    if self.dim_charges is not None:
        for head in self.dim_charges:
            self.dim_charges[head] = self.dim_charges[head][indx]
    if self.dim_pol is not None:
        self.dim_pol = self.dim_pol[indx]
    if self.dim_efficiencies is not None:
        self.dim_efficiencies = self.dim_efficiencies[indx]
    if self.dim_cross_sections is not None:
        self.dim_cross_sections = self.dim_cross_sections[indx]


def collect_timing(self, f, indices):
    '''Collect the timing info.'''
    from datetime import datetime, timedelta

    # Collect the starting time
    if 'RUNTIME' in indices:
        ix = indices['RUNTIME'] + 1
        tp = ' '.join(f[ix].strip().split()[3:]) # The date
        tp += ' ' + f[ix+1].strip().split()[3]     # Add time to date
        self.start = datetime.strptime(str(tp), '%A, %b %d, %Y %H:%M:%S')
    else:
        self._raise_or_pass('Could not find DIM runtime conditions')

def collect_technical(self, f, indices):
    '''Collect technical info, such as where the job was run and termination'''

    # Look above the timing for an error message
    if 'TIMING' in indices:
        ix = indices['TIMING']
        try:
            # Look up to 50 lines before the timing for an error
            ix = next(i for i in range(ix, ix-50, -1) if 'ERROR' in f[i])
        except StopIteration:
            # If no error was found then we terminated correctly
            self.termination = 'NORMAL TERMINATION'
        else:
            # An error was found, save it
            self.termination = f[ix].strip().replace('ERROR: ', '')

    # Find the number of processor
    if 'RUNTIME' in indices:
        ix = indices['RUNTIME']

        # Get the host name
        self.host = f[ix+3].split()[3]
        # Get the number of processors
        self.nprocs = int(f[ix+4].split()[3])

        # Now make the computer names uniform using predefined names
        if self.host is not None:
            for name in self._computer_names:
                if name.lower() in self.host: self.host = name
