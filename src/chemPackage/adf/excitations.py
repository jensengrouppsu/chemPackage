from __future__ import print_function, division
import re
from numpy import array, append, arange
from numpy import zeros, vectorize, argsort
import numpy as np
# deprcatated
# from numpy import fastCopyAndTranspose as fcat
from numpy import where, isnan
from numpy.core.records import fromarrays
from ..f2py import find_equal
from ..constants import EV2HART

def collect_excitations(self, f, indices):
    '''Collect excitations.'''

    try:
        __collect_excitations_block(self, f, indices)
    except IndexError:
        self._raise_or_pass('Error locating excitation blocks')
    else:

        # Find and match the transition dipole moments
        __match_TDM(self, f, indices)

        # Find and match the transitions
        __collect_transitions(self, f, indices)


def collect_excited_state(self, f, indices):
    '''Collect the excited state gradient and the correct excitations block.'''

    # First collect the excited state identifier as '# state'
    for ln in self.key['EXCITEDGO']:
        if 'STATE' in ln.upper():
            sym, num = ln.split()[1:3]
            self.es = ' '.join([num, __symchange(sym)])

    # Now we collect the first excited state gradient, which is based off
    # optimized ground state geometry.  We want this because it can be used
    # for spectroscopic properties.
    self.es_gradient = {}
    if 'EXCITED STATE GRADIENTS' in indices:
        s = indices['EXCITED STATE GRADIENTS']
        # End of the block is a space
        e = next(i for i, x in enumerate(f[s:], s) if not x.strip())
        # Collect the block
        self.es_gradient['total TDDFT gradient'] = (
            array([x.split()[2:5] for x in f[s:e]],dtype=float))
    else:
        self._raise_or_pass('Error locating excited state gradients')

    # Collect the excited state dipole moment, only for an optimization
    if 'GEOMETRY' in self.calctype:
        if 'EXCITED STATE DIPOLE' in indices:
            ix = indices['EXCITED STATE DIPOLE']
            self.es_dipole = array(f[ix].split()[5:8], dtype=float)
        else:
            self._raise_or_pass('Error locating excited state dipole moment')

    # During the excited state geometry optimization process, the excitations
    # are calculated multiple times.  Only collect the first time.
    if 'EXCITATIONS END' in indices:
        e = indices['EXCITATIONS END']
        # Collect the excitations using a limited file
        for n in ('TDM', 'TRANSITIONS', 'EXACT EXCITATIONS',
                  'DAVIDSON EXCITATIONS'):
            try:
                indices[n] = [x for x in indices[n] if x < e]
            except KeyError:
                pass
        try:
            collect_excitations(self, f[:e], indices)
        except IndexError:
            self._raise_or_pass('Error locating excitation blocks')
    else:
        self._raise_or_pass('Error locating the end of the excitation block')


def __collect_excitations_block(self, f, indices):
    '''Collect excitations from the excitations block..'''

    # The excitation energies are printed together in a list, but they
    # are grouped by type.  All singlet-singlet are in a table, singlet-triplet
    # are in a separate table, and spin-unrestricted are in a separate table.
    # All excitations from all tables are grouped together, then sorted by
    # increasing energy.  A separate array keeps track of the excitation type.

    # Find all possible locations of the excitation blocks.
    start = []
    for label in ('SPIN-UNRESTRICTED', 'SINGLET-SINGLET', 'SINGLET-TRIPLET'):
        if label in indices:
            start.append((indices[label], label))

    # If none were found, raise an error
    if not start:
        raise IndexError

    # Map the full label to shorthand
    label_change = { 'SPIN-UNRESTRICTED' : 'UN',
                     'SINGLET-SINGLET'   : 'SS',
                     'SINGLET-TRIPLET'   : 'ST' }

    # Collect all the blocks
    een = []
    osc = []
    sym = []
    lab = []
    for s, label in start:

        # List of excitations ends at first blank
        e = next(i for i, x in enumerate(f[s:], s) if not x.strip())

        # Collect the main excitations block into a numpy array
        for line in f[s:e]:
            ln = line.split()
            # adf2018 changed how the excitation block was printed
            try:
                rando = ln[5]
                try:
                    een.append(float(ln[1]))
                    osc.append(float(ln[3]))
                    sym.append(__symchange(ln[5]))
                    lab.append(label_change[label])
                except ValueError:
                    pass # probably excitation too high in FDE(EO) case

            except IndexError:
                try:
                    een.append(float(ln[1]))
                    osc.append(float(ln[3]))
                    sym.append(__symchange(ln[4]))
                    lab.append(label_change[label])
                except ValueError:
                    pass # probably excitation too high in FDE(EO) case

    # Now, place them in order of increasing energy
    self.nexcite = len(een)
    index = argsort(array(een))
    self.excitation_energies = array(een)[index]
    self.oscillator_strengths = array(osc)[index]
    self.excitation_symmetries = array(sym)[index]
    self.excitation_type = array(lab)[index]

    # Make a private index for excitations
    self._excite_index = arange(self.nexcite, dtype=int)


def __match_TDM(self, f, indices):
    '''Find and match the transition dipole moments to the excitations.'''

    # The transition dipole moments appear in the output file before
    # the excitations.  However, if a transition dipole moment is weak
    # it is not printed at all, so number of transition dipole moments
    # will be less than or equal to the number of excitations.
    # Note that if there was symmetry in the calculation the TDM will be
    # split into separate lists according to symmetry.  Therefore, the
    # index of where each list starts is found and for each list the TDM
    # are collected.  Because the TDM are not necessarily listed in the
    # same order as the excitations the TDM energies and oscillator
    # strengths are compared to the excitations and if both match then
    # the TDM is correlated to that excitation.  If no TDM is listed for
    # an excitation it is defaulted to zero.

    # This may appear more than once if there is symmetry
    if 'TDM' in indices:
        ar = indices['TDM']
    else:
        for subkeys in self.key['EXCITATION']:
            if 'EXACT' in  subkeys.upper():
                return # TDM not printed for exact diagonalization
#        self._raise_or_pass('Error locating transition dipole moment blocks')
        return

    # Find all the transition dipole moments.  There may not be any
    TDM = array([])
    for s in ar:
        e = next(i for i, x in enumerate(f[s:], s) if not x.strip())
        tp = array([ln.split()[1:6] for ln in f[s:e]], dtype=float)
        TDM = append(TDM, tp.flatten())
    TDM = TDM.reshape(-1,5)

    # Default TDM to zero == test_line else False
    self.TDM = zeros((self.nexcite, 3), dtype=float)

    # Match up the TDM to the excitations.
    # If none were found, then all were weak
    try:
        test = TDM[0,0]
        test = TDM[0,1]
    except IndexError:
        pass
    else:
        indx = find_equal(TDM[:,0]*EV2HART, self.excitation_energies,
                          TDM[:,1], self.oscillator_strengths)
        for i, ix in enumerate(indx):
            if ix == -1: continue
            self.TDM[i] = TDM[ix,2:5]


def __collect_transitions(self, f, indices):
    '''Find and collect the restricted transitions, then match them to
    excitations, using a similar algorithm to the TDM.'''

    # The transitions suffer from a similar problem to the TDM, in that
    # they are grouped by symmetry.  There is a 1:1 correlation to the 
    # transitions with the excitation table for each symmetry group, so
    # this table is used to match the transitions to the excitations in
    # the same way as the TDM (i.e. using find_equal).

    # Note that each transition table is space-separated.

    # Vectorize the transition symmetry conversion function to be broadcastable
    def tsym(x):
        '''Function to format trans sym groups.'''
        m = re.search(r'(\d+)(.*)', x)
        return m.group(1)+' '+__symchange(m.group(2).capitalize())
    trans_sym = vectorize(tsym)

    def tbend(x, y, z):
        '''Function to determine the end of a transition block.'''

        test_line = 'Eigenvalues of small (approximate) problem'
        # If the current line is blank, then it may be the end
        if not x.strip():
            # If the next line is also blank or next line has the test line,
            # then this is the end. Otherwise continue.
            #return True if (y == test_line) or (len(z.split())<7 and z!='') else False
            # Zhongwei: make it work for 200 lowest excited states of au25
            return True if (y == test_line) or (len(z.split())<5 and z!='') else False

        # If the current line is not blank, then continue
        else:
            return False
        
    # First find where the excitation lists above the transitions are
    # Different between excact and Davidson
    if 'EXACT' in self.subkey:
        if 'EXACT EXCITATIONS' in indices:
            ar = indices['EXACT EXCITATIONS']
        else:
            self._raise_or_pass('Error locating exact excitation energies')
            return
    else:
        if 'DAVIDSON EXCITATIONS' in indices:
            ar = indices['DAVIDSON EXCITATIONS']
        else:
            self._raise_or_pass('Error locating davidson excitation energies')
            return
    # Next find where the lists of transitions are
    # Place the locations of these two things into one array together
    if 'TRANSITIONS' in indices:
        # ar = fcat(array([ar,indices['TRANSITIONS']]))
        ar = np.copy(np.array([ar, indices['TRANSITIONS']]).T)
    else:
        self._raise_or_pass('Error locating electronic transitions')
        return

    # Loop over each transition block
    trans = []
    en  = array([], dtype=float)
    osc = array([], dtype=float)

    for ix in ar:

        # Collect osc strength and energy that corresponds to each transition
        s = ix[0]
        e = next(i for i,x in enumerate(f[s:], s) if not x.strip())
        tp = array([ln.split() for ln in f[s:e]])
        en = append(en, array(tp[:,1], dtype=float))
        osc = append(osc, array(tp[:,3], dtype=float))

        # Collect the transitions themselves
        # First Line of this transition block
        s = ix[1]
        s = next(i for i, x in enumerate(f[s:], s) if not x.strip())
        if not f[s+1].strip(): s += 1
        # End if this block this matches the break rules 
        e = next(i for i, x in enumerate(f[s:], s) if tbend(x, f[i+1], f[i+2])) + 1

        # Now find where all the blanks are in the current block
        iend = [i for i, x in enumerate(f[s:e], s) if x == '']

        # Remove double blanks for newer version of ADF
        doubles = []
        for i in range(len(iend)-1):
            if iend[i] == iend[i+1]-1: doubles.append(iend[i])
        for i in doubles:
            iend.remove(i)

        # Collect from blank to blank
        for i in range(0, len(iend)-1, 1):

            # Sometimes, if the transition is weak enough, the strength is left
            # blank.  Account for this by adding zeros.  Then collect
            ar = []
            for ln in f[iend[i]+1:iend[i+1]]:
                if ln == '': continue
                tp = ln.split()
                # Should be 8 (9) elements... if 5 (6), fill to 8 (9)
                if 'UNRESTRICTED' in self.calctype:
                    if len(tp) == 6: tp.extend([0.0, 0.0, 0.0])
                else:
                    if len(tp) == 5: tp.extend([0.0, 0.0, 0.0])
                ar.append(tp)

            # Store in the correct locations
            ar = array(ar)
            add = 1 if 'UNRESTRICTED' in self.calctype else 0
            occ = trans_sym(ar[:,1+add])
            unocc = trans_sym(ar[:,3+add])
            pcent = array(ar[:,4+add], dtype=float)
            # Spin is 1 for alpha, 2 for beta, or 0 for restricted
            if 'UNRESTRICTED' in self.calctype:
                spin = array([{ 'Alph' : 1, 'Beta' : 2 }[x] for x in ar[:,1]])
            else:
                spin = zeros(len(ar[:]), dtype=int)
            # Add trans as a record array - numpy equivalent to dictionary
            trans.append(fromarrays([occ, unocc, pcent, spin],
                                    names='occ,unocc,pcent,spin'))

    # Now match the transitions to the excitations
    indx = find_equal(en, self.excitation_energies,
                      osc, self.oscillator_strengths)

    try:
        self.transitions = []
        for ix in indx:
            if ix == -1: continue
            self.transitions.append(trans[ix])
        self.transitions = tuple(self.transitions)
    except IndexError:
        pass


def __symchange(sym):
    '''Function to standardize the summetry group.'''
    m = re.match(r"([A-Z])(\d?)('*)$", sym)
    try: # Convert i.e. E2'' => EEE2
        return m.group(1)+m.group(1)*len(m.group(3))+m.group(2)
    except (IndexError, AttributeError):
        return sym

def collect_fdec_excitations(self, f, indices):
    '''Collects the excitations and transition dipoles from an FDEc
    subsystem excitation calculation.'''

    if 'FDEC ENERGIES' in indices:
        tempe = array([], dtype=float)
        tempt = array([], dtype=float)
        ar = indices['FDEC ENERGIES']
        for ix in ar:
            for i in range(ix,len(f)):
                line = f[i]
                if line=='': break
                if len(line.split()) == 3:
                    tempe = append(tempe, line.split()[1])
                    tempt = append(tempt, line.split()[2])
        self.excitation_energies = EV2HART(array(tempe, dtype=float))
        self.oscillator_strengths = array(tempt, dtype=float)

        # change NaN to 0.
        indx = where(isnan(self.excitation_energies))
        self.excitation_energies[indx] = 0.
        indx = where(isnan(self.oscillator_strengths))
        self.oscillator_strengths[indx] = 0.
