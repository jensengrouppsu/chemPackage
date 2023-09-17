from __future__ import print_function, division
from numpy import array, where, vectorize, zeros, argsort, append
from numpy.core.records import fromarrays
import re

def collect_orbitals(self, f, indices):
    '''Collect orbital information from the orbital block.'''

    # The orbital information is given in a list from lowest energy MO
    # to highest energy MO.  Within that list, it gives the AOs that
    # make up the MO in order of percent contribution.  Each line with
    # MO information is outdented, and the AO with the largest
    # contribution is on the same line.
    # For unrestricted calculations, there are two such lists.

    # Symmetry group transformation
    transform = {'A1.g': 'S+.g', 'A2.u': 'S+.u', 'E1.g': 'Pi.g',
                 'E1.u': 'Pi.u', 'E2.g': 'De.g', 'E2.u': 'De.u'}

    # Start at the beginning of the orbital blocks
    start = indices['ORBITALS']

    # Initiallize
    aos = []
    self.orbital_energies = array([], dtype=float)
    self.orbital_occ = array([], dtype=float)
    self.orbital_ids = array([], dtype=str)
    self.orbital_spin = array([], dtype=int)

    # loop over each orbital block
    for n, s in enumerate(start):

        # Location of end of orbital block.  It is a blank line, a line with
        # only the number one on it, or conatins 'Using NEW gradient routines'
        fn = lambda x: True if not x.strip(' 1') or 'NEW' in x or len(x) > 95 else False
        e = next(i for i, x in enumerate(f[s:], s) if fn(x))

        # Exit now if there are no orbitals
        if s == e: return

        # Loop over the orbital block, and find all the lines that have MOs
        # These are the lines with 11 fields.
        imo = [i for i, x in enumerate(f[s:e], s) if len(x.split()) == 11 or len(x.split()) == 12]

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
    index = argsort(self.orbital_energies)
    self.orbital_energies = self.orbital_energies[index]
    self.orbital_occ = self.orbital_occ[index]
    self.orbital_ids = self.orbital_ids[index]
    try:
        self.orbital_spin = self.orbital_spin[index]
    except IndexError:
        pass
    self.atomic_orbitals = []
    for ix in index:
        self.atomic_orbitals.append(aos[ix])

    # Convert to tuple
    self.atomic_orbitals = tuple(self.atomic_orbitals)

    # Determine HOMO-LUMO info
    __find_homo_lumo(self)


def __collect_mos(self, f, imo):
    '''Collect the MOs in the MO block.'''

    # For each MO, the energy is the first field, the self.orbital_occupation is the second,
    # and the ID is the third and fourth.

    # Symmetry group transformation
    transform = {'A1.g': 'S+.g', 'A2.u': 'S+.u', 'E1.g': 'Pi.g',
                 'E1.u': 'Pi.u', 'E2.g': 'De.g', 'E2.u': 'De.u'}
    # Degeneracy remover
    deg = re.compile(r'(.+?):\d')

    # Allocate
    orb_en  = []
    orb_occ = []
    orb_ids = []

    # Loop over each line with the MO
    for ln in [f[x] for x in imo]:
        l = ln.split()
        orb_en.append(float(l[0]))
        orb_occ.append(float(l[1]))
        # Remove degeneracy from ID, then transform if necessary
        mo_id = deg.sub(r'\1', l[3])
        orb_ids.append(l[2]+' '+transform.get(mo_id, mo_id))

    # Convert to numpy arrays and return
    return array(orb_en), array(orb_occ), array(orb_ids)


def __collect_aos(self, f, imo, e):
    '''Collect the AOs in the MO block.'''

    # Each AO list for an MO starts in the same line as an MO, then
    # continues to the line before the next MO.

    # Function to split lines containing AO info and return
    def ao_lines(line):
        ln = line.split()
        try:
            pcent = float(ln[0].rstrip('%')) / 100
        except ValueError:
            pcent = float('nan')
        try:
            ao_id = ln[5]+' '+ln[6]
            sym = ln[1]+' '+ln[2]
        except IndexError:
            if ln[1].upper() == 'CORE':
                ao_id = 'core'
                sym = 'core'
            else:
                raise IndexError
        return pcent, ao_id, sym
    # Vectorize the above function to work on numpy arrays as a whole
    aolines = vectorize(ao_lines)

    # Loop over each MO index location, which is also the start of the AO list
    aos = []
    #for i in xrange(len(imo)):
    for i in range(len(imo)):
        # Find the index range for each set of AOs.
        # The try is if this is the last MO then imo[i+1] does not exist
        try:
            irange = range(imo[i], imo[i+1])
        except IndexError:
            irange = range(imo[i], e)
        lines = [f[x] for x in irange]
        # The first index needs to be modified to remove the MO information.
        lines[0] = ' '.join(lines[0].split()[4:])
        # Now split and collect the info from each line,
        try:
            pcent, ao_id, sym = aolines(lines)
            # Place into a record array
            aos.append(fromarrays([pcent, ao_id, sym], names='pcent,ao_id,sym'))
        except IndexError:
            pass

    return tuple(aos)


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
