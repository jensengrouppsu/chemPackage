from __future__ import print_function, division
from numpy import array, where, zeros, argsort, append
from numpy.core.records import fromarrays
from prep import index_natsorted
import re

def collect_orbitals(self, f, indices):
    '''Collect orbital information from the orbital block.'''

    # The orbital information is given in large block broken into several
    # smaller sub-blocks.  Each of these sub-blocks begins with the
    # information on the MO, then contains a two-column list of the
    # AOs in that MO.

    try:
        __collect_orbital_block(self, f, indices)
    except IndexError:
        self._raise_or_pass('Error collecting MO analysis')
    else:
        # Determine HOMO-LUMO info
        __find_homo_lumo(self)


def __collect_orbital_block(self, f, indices):
    '''Collect orbital information from the orbital block.'''

    if 'RESTRICTED DFT' in indices:
        start = [indices['RESTRICTED DFT']+3]
    elif 'UNRESTRICTED DFT ALPHA' in indices:
        try:
            start = [indices['UNRESTRICTED DFT ALPHA']+3,
                     indices['UNRESTRICTED DFT BETA']+3]
        except KeyError:
            self._raise_or_pass('Error locating DFT BETA orbitals')
    elif 'RESTRICTED HF' in indices:
        start = [indices['RESTRICTED HF']+3]
    elif 'UNRESTRICTED HF ALPHA' in indices:
        try:
            start = [indices['UNRESTRICTED HF ALPHA']+3,
                     indices['UNRESTRICTED HF BETA']+3]
        except KeyError:
            self._raise_or_pass('Error locating UHF orbitals')
    else:
        self._raise_or_pass('Error locating ground state molecular orbitals')

    # Allocate
    self.orbital_energies = array([], dtype=float)
    self.orbital_occ = array([], dtype=float)
    self.orbital_ids = array([], dtype=str)
    self.orbital_spin = array([], dtype=int)
    aos = []

    # Loop over each spin
    for n, s in enumerate(start):

        # Find end of orbital block
        # Newer versions of NWChem print alpha-beta electron overlaps.
        for i, x in enumerate(f[s:], s):
            if 'alpha - beta orbital overlaps' in x:
                e = i - 3
                break
            if 'center of mass' in x:
                e = i - 2
                break
            if 'Total Density - Mulliken Population Analysis' in x:
                e = i - 3
                break
            if 'DFT Final Beta Molecular Orbital Analysis' in x:
                e = i - 2
                break
            if 'UHF Final Beta Molecular Orbital Analysis' in x:
                e = i - 2
                break

        # Locate each MO in the block
        imo = [i for i, x in enumerate(f[s:e], s) if 'Vector' in x]

        # Collect the MOs
        mo_info = __collect_mos(self, f, imo)
        # Collect the AOs
        aos.extend(__collect_aos(self, f, imo, e))

        # Assign the MO info
        self.orbital_energies = append(self.orbital_energies, mo_info[0])
        self.orbital_occ = append(self.orbital_occ, mo_info[1])
        self.orbital_ids = append(self.orbital_ids, mo_info[2])
        # 1 for spin up, 2 for spin down, 0 for restricted
        if 'UNRESTRICTED' in self.calctype:
            spin = zeros(len(mo_info[0]), dtype=int) + ( n + 1 )
            self.orbital_spin = append(self.orbital_spin, spin)
        else:
            self.orbital_spin = zeros(len(mo_info[0]), dtype=int)

    # Sort by increasing energy.  Necessary because of unrestricted possibility
    self.nmos = len(self.orbital_energies)
    index = argsort(self.orbital_energies, kind='mergesort')
    index2 = index_natsorted(self.orbital_ids[index])
    self.orbital_energies = self.orbital_energies[index][index2]
    self.orbital_occ = self.orbital_occ[index][index2]
    self.orbital_ids = self.orbital_ids[index][index2]
    self.orbital_spin = self.orbital_spin[index][index2]
    self.atomic_orbitals = []
    for ix in index[index2]:
        self.atomic_orbitals.append(aos[ix].copy())
    self.atomic_orbitals = tuple(self.atomic_orbitals)

def __collect_mos(self, f, imo):
    '''Collect the MOs in the MO block.'''

    # Allocate
    orb_en  = []
    orb_occ = []
    orb_ids = []

    # Loop over each line with the MOs
    for ln in [f[i] for i in imo]:
        # The occupation number. Round to two decimals.
        # It is print off as 1D3 insead of 1E3, so change D to E
        orb_occ.append(round(float(ln[18:30].replace('D', 'E')), 2))
        # The MO energy
        orb_en.append(float(ln[34:47].strip().replace('D', 'E')))
        # ID is # plus symmetry.  If symmetry is not printed, it is A
        try:
            sym = ln[58:].strip().upper()
        except IndexError:
            sym = 'A'
        # Attach the symmetry group to the orbital number
        orb_ids.append(' '.join([ln[7:12].strip(), sym]))

    # Convert to numpy arrays and return
    return array(orb_en), array(orb_occ), array(orb_ids)


def __collect_aos(self, f, imo, end):
    '''Collect the AOs in the MO block.'''
    from mfunc import norm, mult
    from numpy import dot

    # Initiallize the AO list
    aos = []

    # Loop over each MO sub-block
    for i, start in enumerate(imo):
        
        # Initiallize arrays for this MO
        pcent = array([], dtype=float)
        ao_id = array([], dtype=str)
        sym = array([], dtype=str)
        cont_sum = 0.0

        # Define start
        s = start + 4

        # Define end
        try:
            e = imo[i+1] - 1
        except IndexError:
            e = end

        # Loop over the AO lines for this MO
        for x in f[s:e]:
            # Atomic orbital number and symmetry
            s = x[22:36].split()
            ao_s = s[0] + ' ' + s[1]
            ao_id = append(ao_id, ao_s)
            # Coefficient in MO
            coeff = float(x[8:19].strip())
            pcent = append(pcent, coeff)
            # The Symmetry
            sym_s = ''.join(s[2:])
            sym = append(sym, sym_s.strip())
            # If two columns then lather, rince and repeat
            if len(x.strip()) > 38:
                s = x[63:77].split()
                ao_s = s[0] + ' ' + s[1]
                ao_id = append(ao_id, ao_s)
                coeff = float(x[49:60].strip())
                pcent = append(pcent, coeff)
                sym_s = ''.join(s[2:])
                sym = append(sym, sym_s)

        # Normalize
        pcent = (pcent/norm(pcent))**2

        # Place into the ao list
        aos.append(fromarrays([pcent, ao_id, sym], names='pcent,ao_id,sym'))

    return aos


def __find_homo_lumo(self):
    '''Determine the HOMO and LUMO based on self.orbital_occupation numbers.'''
    # LUMO is first obital where self.orbital_occupation equals 0.
    # HOMO is last non-zero orbital
    # In cases where there are no virtual orbitals, there is no LUMO
    oe = self.orbital_energies
    oi = self.orbital_ids
    oo = self.orbital_occ
    os = self.orbital_spin
    if 'UNRESTRICTED' in self.calctype:
        self.HOMO = []
        self.LUMO = []
        for i in (0, 1):
            ispin = where(os == i+1)[0]
            try:
                ixl = where(oo[ispin] == 0.0)[0][0]
            except IndexError:
                self.LUMO.append(('', 0.0))
                ixl = self.nmos
            else:
                self.LUMO.append((oi[ispin][ixl], oe[ispin][ixl], ixl))
            ixh = where(oo[ispin] > 0.0)[0][-1]
            self.HOMO.append((oi[ispin][ixh], oe[ispin][ixh], ixh))

            # Issue a HOMO - LUMO warning
            if ixl < ixh:
                self._raise_or_pass('LUMO energy less than the HOMO')
            elif ixl != ixh + 1:
                self._raise_or_pass('HOMO and LUMO not sequential')
        # Make into a tuple
        self.HOMO = tuple(self.HOMO)
        self.LUMO = tuple(self.LUMO)

    else:
        try:
            ixl = where(oo == 0.0)[0][0]
        except IndexError:
            self.LUMO = ('', 0.0)
            ixl = self.nmos
        else:
            self.LUMO = (oi[ixl], oe[ixl], ixl)
        ixh = where(oo > 0.0)[0][-1]
        # If there was a 1.0 occ, this was a SOMO
        if oo[ixh] == 1.0:
            self.SOMO = (oi[ixh], oe[ixh], ixh)
            self.HOMO = (oi[ixh-1], oe[ixh-1], ixh-1)
        else:
            self.HOMO = (oi[ixh], oe[ixh], ixh)

        if ixl < ixh:
            self._raise_or_pass('LUMO energy less than the HOMO')
        elif ixl != ixh + 1:
            self._raise_or_pass('HOMO and LUMO not sequential')
