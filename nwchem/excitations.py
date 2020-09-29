from __future__ import print_function, division
from ..constants import EV2HART
from numpy import array, where, argsort, arange
from numpy.core.records import fromarrays
import numpy

def collect_excitations(self, f, indices):
    '''Driver to collect excitations.'''
    # There are different headers for different levels of theory
    if 'DFT EXCITATIONS' in indices:
        __dft_excitations(self, f, indices['DFT EXCITATIONS'])
    elif 'CCSD EXCITATIONS' in indices:
        __ccsd_excitations(self, f, indices['CCSD EXCITATIONS'])
    # Add more drivers here for other methods
    #elif 'OTHER METHOD' in indices:
    else:
        self._raise_or_pass('No excitation block was found')


def __dft_excitations(self, f, ix):
    '''Collect the DFT excitations'''

    # Map the full label to shorthand
    label_change = { 'singlet' : 'SS',
                     ' single' : 'SS',
                     'doublet' : 'SD',
                     'triplet' : 'ST' }

    # Table for converting spin-multiplicity to a python index.
    mult_map = { '1' : 'singlet',
                 '2' : 'doublet',
                 '3' : 'triplet',}

    # Initiallize lists
    een = []
    osc = []
    sym = []
    lab = []
    TDM = []
    TQM = []
    trans = []
    if "CD SPECTRUM" in self.calctype:
        MDM = [] # Magnetic transition dipole
        ors = [] # Optical rotatory strength

    # Locate where each excitation is given
    arr = [i for i, x in enumerate(f[ix+1:], ix+1) 
                if 'Root' in x and 'eV' in x]

    # Loop over each excitation
    for s in arr:

        # Collect the information on the Exciation line
        # Use the strict format of NWChem to do column-based splitting.
        # singlet/triplet, symmetry, energy
        #ln = (f[s][10:17], f[s][17:24], f[s][24:35])
        
        # Using python split function now since column-based splitting
        # is unreliable when NWChem changes its formatting from version
        # to version
        ln = f[s].split()
        #print (ln)
        # Energy
        if len(ln) == 8:
            een.append(float(ln[4]))
        elif len(ln) == 7:
            een.append(float(ln[3]))
        elif len(ln) == 9:
            een.append(float(ln[4]))
        # Type (singlet/triplet).  
        if 'RESTRICTED' in self.calctype:
            lab.append(label_change[ln[2]])
        else:
            mult = mult_map[str(self.spin_multiplicity)]
            lab.append(label_change[mult])
        # Symmetry.
        if len(ln) == 8:
            sym.append(ln[3].strip().upper())
        elif len(ln) == 7:
            sym.append(ln[2].strip().upper())
        elif len(ln) == 9:
            sym.append(ln[3].strip().upper())

        # Collect TDMs.  Try is in case of forbidden excitations
        if 'RESTRICTED' in self.calctype:
            # Restricted calculation
            ln = f[s+2].split()
        else:
            # Unrestricted calculation
            ln = f[s+3].split()
        try:
            TDM.append([float(ln[3]), float(ln[5]), float(ln[7])])
        except IndexError:
            TDM.append([0.0, 0.0, 0.0])
        except ValueError:
            # This is needed for triplet calculations
            TDM.append([0.0, 0.0, 0.0])

        # Collect oscillator strengths.  Newer versions of NWChem
        # place this differently. 

#       This first line used to work, but causes issues with NWChem dev
        if 'Dipole Oscillator Strength' in f[s+5]:
            osc.append(float(f[s+5].split()[3]))
            if 'Occ.' in f[s+7]:
                st = s + 7 # Where the transitions start
            else:   # This was needed for nwchem-dev 07-2016
                st = s + 10
        elif 'Dipole Oscillator Strength' in f[s+9]:
            osc.append(float(f[s+9].split()[3]))
            st = s + 11 # Where the transitions start
        elif 'Dipole Oscillator Strength' in f[s+6]:
            # For unrestricted calculations
            osc.append(float(f[s+6].split()[3]))
            st = s + 8
        elif 'forbidden' in f[s+3]:
            # For triplet calculations
            osc.append(float(0.0))
            st = s + 5
        elif 'Oscillator Strength' in f[s+3]:
            # For really old restricted NWChem calculations
            osc.append(float(f[s+3].split()[2]))
            st = s + 5
        elif 'Oscillator Strength' in f[s+4]:
            # For really old unrestricted NWChem calculations
            osc.append(float(f[s+4].split()[2]))
            st = s + 6
        else:
            # Try is in case of forbidden excitations
            try:
                osc.append(float(f[s+3].split()[2]))
            except ValueError:
                osc.append(float(0.0))
            st = s + 5

        # Collect quadrupole transition moments
        if 'RESTRICTED' in self.calctype:
            ln1 = f[s+3][:]
            ln2 = f[s+4][:]
        else:
            ln1 = f[s+3][:]
            ln2 = f[s+4][:]
        try:
            TQM.append([float(ln1[28:38]), float(ln1[41:51]), float(ln1[54:64]),
                        float(ln1[28:38]), float(ln1[41:51]), float(ln1[54:64])])
        except ValueError:
            TQM.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        except IndexError:
            TQM.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

        # Collect magnetic transition dipole moments and optical rotatory 
        # strengths for Circular Dichroism.
        if 'CD SPECTRUM' in self.calctype:
            if 'RESTRICTED' in self.calctype:
                ln = f[s+10].split()
            else:
                ln = f[s+11].split()
            MDM.append([float(ln[1]),float(ln[3]),float(ln[5])])
            # Here, it is assumed the velocity gauge is more accurate
            # than the length gauge for optical activity
            if 'RESTRICTED' in self.calctype:
                if 'VELOCITY' in self.calctype:
                    ors.append(float(f[s+17].split()[5]))
                else:
                    ors.append(float(f[s+11].split()[4]))
            else:
                if 'VELOCITY' in self.calctype:
                    ors.append(float(f[s+18].split()[5]))
                else:
                    ors.append(float(f[s+12].split()[4]))

        # If Circular Dichroism is calculated, the transitions start
        # later in the output file.
        if 'VELOCITY' in self.calctype:
            if 'Rotatory Strength' in f[s+17]:
                st = s + 19 # For restricted calculations
            elif 'Rotatory Strength' in f[s+18]:
                st = s + 20 # For unrestricted calculations
        else:
            if 'Electric Transition Dipole' in f[s+7]:
                st = s + 13 # For restricted calculations
            elif 'Electric Transition Dipole' in f[s+8]: 
                st = s + 14 # For unrestricted calculations

        # Find where the transitions end.  Blank or series of dashes.
        fn = lambda x: True if not x.strip() or '----' in x else False
        e = next(i for i , x in enumerate(f[st:], st) if fn(x))

        # Collect transitions.  Use the strict format of NWChem
        # to do column-based splitting.  Order is:
        # occupied start, occupied end
        # unoccupied start, unoccupied end
        # percent
        # Updated to use python split function since formatting changes
        # from version to version with nwchem.
        try:
            # New versions of NWChem are formatted as Occ. ###
            if 'RESTRICTED' in self.calctype:
                try:
                    ar = array([[x.split()[1], x.split()[2].upper(),
                                 x.split()[5], x.split()[6].upper(),
                                 x.split()[7]] for x in f[st:e]])
                except IndexError:
                    ar = array([[x[9:12], x[14:16].upper(),
                                 x[28:31], x[32:34].upper(),
                                 x[36:46]] for x in f[st:e]])
            else:
                try:
                    ar = array([[x.split()[1], x.split()[3].upper(),
                                 x.split()[6], x.split()[8].upper(),
                                 x.split()[9]] for x in f[st:e]])
                except IndexError:
                    ar = array([[x[9:12].strip(), x[17:24].strip().upper(),
                                 x[33:36].strip(), x[41:46].strip().upper(),
                                 x[47:56].strip()] for x in f[st:e]])
            
        except(ValueError):
            #Old versions of NWChem were formatted as Occ.###
            if 'RESTRICTED' in self.calctype:
                ar = array([[x[9:12].strip(), x[12:18].strip().upper(),
                            x[28:31].strip(), x[31:36].strip().upper(),
                            x[36:45].strip()] for x in f[st:e]])
            else:
                ar = array([[x[9:12].strip(), x[17:24].strip().upper(),
                            x[33:36].strip(), x[41:46].strip().upper(),
                            x[47:56].strip()] for x in f[st:e]])
        # Determine occupied then unoccupied orbitals in transitions
        # Number, then symmetry group
        occ   = array([' '.join([x, y]) for x, y in zip(ar[:,0], ar[:,1])])
        unocc = array([' '.join([x, y]) for x, y in zip(ar[:,2], ar[:,3])])
        # Determine the percent contribution.
        pcent = array(ar[:,4], dtype=float)**2
        # Normalize to the sum of squares
        pcent = (pcent/numpy.linalg.norm(pcent))**2
        # Sort according to decreasing percent
        index = argsort(pcent)[::-1]
        # Add to the transitions list as a record array
        trans.append(fromarrays([occ[index], unocc[index], pcent[index]],
                                names='occ,unocc,pcent'))
    # Now that we've collected everything, sort the energies.
    # This is required in case we calculated both singlets and triplets
    self.nexcite = len(een)
    index = argsort(array(een))
    self.excitation_energies = array(een)[index]
    self.oscillator_strengths = array(osc)[index]
    self.excitation_symmetries = array(sym)[index]
    self.excitation_type = array(lab)[index]
    self.TDM = array(TDM)[index]
    self.TQM = array(TQM)[index]
    self.transitions = []
    for ix in index:
        self.transitions.append(trans[ix])
    self.transitions = tuple(self.transitions)

    # Circular Dichroism related stuff
    if 'CD SPECTRUM' in self.calctype:
        self.MDM = array(MDM)[index]
        self.opt_rot_strengths = array(ors)[index]

    # Make a private index for excitations
    self._excite_index = arange(self.nexcite, dtype=int)

def __ccsd_excitations(self, f, ix):
    '''Collect the CCSD excitations.'''

    # Initiallize lists
    een = []
    osc = []
    sym = []
    lab = []
    TDM = []
    trans = []

    # Locate each excitation
    arr = [i for i, x in enumerate(f[ix-4:], ix-4) if 'Excited state' in x]

    # Loop over each excitation
    for s in arr:
        # Energy
        een.append(float(f[s+1].split('=')[-1].strip()))
        # Collect TDMs
        ln = f[s+6].split()
        TDM.append([float(ln[1]), float(ln[3]), float(ln[5])])
        # Collect oscillator strengths
        osc.append(float(f[s+7].split()[-1]))
        # Type (not sure what is printed for triplets)
        lab.append('SS')
        # Symmetry (doesn't seem to print anything for symmetry)
        sym.append('A')

    # Now that we've collected everything, sort the energies.
    # This is required in case we calculated both singlets and triplets
    self.nexcite = len(een)
    index = argsort(array(een))
    self.excitation_energies = array(een)[index]
    self.oscillator_strengths = array(osc)[index]
    self.excitation_symmetries = array(sym)[index]
    self.excitation_type = array(lab)[index]
    self.TDM = array(TDM)[index]

    # Make a private index for excitations
    self._excite_index = arange(self.nexcite, dtype=int)

def collect_excited_state_gradient(self, f, indices):
    '''Collect the excited state gradient for a TDDFT calculation.'''

#    # First collect the excited state identifier 
#    for ln in self.key['GRAD']:
#        # This is coded in a silly way, I need to fix this.  We assume
#        # no symmetry for the moment.  Eventually, the gradients may
#        # be written in a more clever way.
#        if 'ROOT' in ln.upper():
#            num = ln.split()[1]
#            self.es = ' '.join([num, 'A'])

    # Initialize the gradient dictionary and location dictionary
    self.es_gradient = {}

    # Find the gradients
    if 'TDDFT GRAD START' in indices:
        s = indices['TDDFT GRAD START'][-1]

        # Collect the excited state identifier.
        # NB, because of how it's currently printed,
        # we have to assume no symmetry.
        for ix in range(s-5,s):
            if 'ROOT' in f[ix].upper():
                num = f[ix].split()[1]
                self.es = ' '.join([num, 'A'])

        if 'TDDFT GRAD END' in indices:
            e = indices['TDDFT GRAD END'][-1]
        else:
            e = indices['FILE END']
        arr = [i for i, x in enumerate(f[s:e], s+1)
                    if 'TDDFT ENERGY GRADIENTS' in x]
        for ix in arr:
            self.es_gradient['total TDDFT gradient'] = ( 
                array([x.split()[5:8] for x in f[ix+3:ix+3+self.natoms]], dtype=float))
    else:
        self._raise_or_pass('Could not find excited state gradient')    

    # Old code
    ## Initialize the gradient dictionary and location dictionary
    #self.es_gradient = {}
    #esg_locations = {}

    ## Table of things to search for
    #table = ('Vxc gradient',
    #         'CD gradient',
    #         'fxc gradient',
    #         'nuclear repulsion gradient',
    #         'weighted density gradient',
    #         'kinetic energy gradient',
    #         '2-electron gradient',
    #         'TDDFT CD+XC gradient',)

    ### Find the gradients
    #if 'TDDFT GRAD START' in indices:
    #    s = indices['TDDFT GRAD START'][-1]
    #    e = indices['DFT GRAD END'][-1]
    #    count = s
    #    for ln in f[s:e]:
    #        # Find locations of contributions to the gradient (if they're printed)
    #        for item in table:
    #            if item in ln:
    #                esg_locations[item] = count
    #        # Collect the TDDFT excitation energy gradient
    #        if 'TDDFT EXCITATION ENERGY GRADIENTS' in ln:
    #            self.es_gradient['TDDFT excitation energy gradient'] = ( 
    #                array([x.split()[5:8] for x in f[count+4:count+4+self.natoms]], dtype=float))
    #        # Collect the TDDFT excited state gradient
    #        if 'TDDFT EXCITED STATE ENERGY GRADIENTS' in ln:
    #            self.es_gradient['total TDDFT gradient'] = (
    #                array([x.split()[5:8] for x in f[count+4:count+4+self.natoms]], dtype=float))
    #        count += 1
    #    # Collect the contributions to the gradient
    #    for item in esg_locations.keys():
    #        s = esg_locations[item] + 1
    #        e = esg_locations[item] + 1 + self.natoms
    #        self.es_gradient[item] = array([x.split()[0:3] for x in f[s:e]], dtype=float)
    #else:
    #    self._raise_or_pass('Could not find excited state gradient.')
