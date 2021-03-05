from __future__ import print_function, division
import re
from numpy import array, append, where
from ..chemdata import ChemData

# Create the NWChem class to grab and hold all the information in the file
class NWChem(ChemData):
    '''Read and conatin data from NWChem files
    
    This class is used to read all information from either an NWChem input
    or output file.  It stores this data in formats suitable for processing.  

    The only required argument for collection is the name of the file.
 
    The data available is listed below.  If the data did not exist in the
    NWChem file, then the variable will default to None (or an empty set
    for subkey and calctype, or empty dict for key).
    
    '''
    #__doc__ = __doc__ + attribute_explanations.format('NWChem')

    # Instantiate the NWChem class
    def __init__(self, filename):
        '''Set up the NWChem class. The only argument is the filename'''
        import os

        ChemData.__init__(self)
        self.program = 'NWChem'
        # Find extention
        ftype = os.path.splitext(filename)[1]
        if ftype not in ('.out', '.nw'):
            raise ValueError (ftype+' not a recognized extention')
        self.filetype = ftype[1:]
        self.filename = filename
        # Several NWChem Keywords. Just a reference.
        # BASIS, CCSD, CHARGE, DFT, DPLOT, ECCE_PRINT, ECHO, ECP, ESP,
        # GEOMETRY, GRADIENTS, MEMORY, MP2 MCSCF, NWPW, PRINT, PROPERTY,
        # PYTHON, RELATIVISTIC, SCF, SCRATCH_DIR, SET, SO, STOP, TASK, TCE,
        # TDDFT, TITLE, UNSET, VSCF

        # Create tuples of general keywords
        self.mixedkeys = ('BASIS', 'GEOMETRY',)
        self.blockkeys = ('CCSD', 'COSMO', 'DFT', 'DPLOT', 'ECP', 'ESP',
                          'FREQ', 'FREQUENCY', 'GRADIENTS', 'HESSIAN', 'MP2',
                          'MCSCF', 'NWPW', 'PROPERTY', 'PYTHON',
                          'RELATIVISTIC', 'SCF', 'SO', 'TCE', 'TDDFT', 'VIB',
                          'VSCF','DIMQM', 'DIMPAR', 'GRAD') 
        #self.linekeys = ('CHARGE', 'ECCE_PRINT', 'MEMORY', 'MULT', 'PRINT',
#   jbecca: remove mult single line option, it should be in the dft block
        self.linekeys = ('CHARGE', 'ECCE_PRINT', 'MEMORY', 'PRINT',
                         'SCRATCH_DIR', 'SET', 'START', 'TASK', 'TITLE',
                         'UNSET','EFIELD')
 #       self.singlekeys = ('ECHO', 'ODFT', 'STOP')
        self.singlekeys = ('ECHO', 'STOP')


    def _collect(self, abort=False):
        '''Collect the NWChem data from file.  *abort* will cause collection
        to stop on an error, instead of ignoring the error.
        
        '''

        # Save the abort status
        self._abort = abort

        # Conventions used in _collect and helper functions:
        # fi = input block
        # fo = output information block
        # ft = timing information block
        # fp = parallel information block
        # s = start of current collection block
        # e = end of current collection block
        # tp = general temporary storage
        # tl = temporary line storage
        # ln = a split line
        # ix = temporary index array
        # ar = temporary array
        # bl = a temporary boolean
        # sl = a line to search for
        # fn = a temporary function
        # m = the results of a regex search

        # Read in file
        from .read_file import read_file
        f, indices = read_file(self)

        # Collect input block
        from .input_block import collect_input
        collect_input(self, f, indices)
        # Determine calculation type
        self.__det_calc_type()

        # Exit now if am input file
        if self.filetype == 'nw': return

        # Collect technical information
        from .technical import collect_technical
        collect_technical(self, f, indices)
        
        self.__collect_energy(f, indices)

        # Collect the dipole moment
        self.__collect_dipole(f, indices)

        # Collect symmetry if it exists
        self.__collect_symmetry(f, indices)

        # Collect charges if they exist
        self.__collect_charges(f, indices)

        # Collect coordinates
        from .coordinates import collect_geometry
        collect_geometry(self, f, indices)
        if 'GEOMETRY' in self.calctype:
            from .coordinates import collect_optimized_geometry
            collect_optimized_geometry(self, f, indices)

        # Collect ground state gradients
        if 'GEOMETRY' in self.calctype and 'EXCITED STATE' not in self.calctype:
            from .coordinates import collect_ground_state_gradient
            collect_ground_state_gradient(self, f, indices)

        # Collect excited state gradients
        if 'GEOMETRY' in self.calctype and 'EXCITED STATE' in self.calctype:
            from .excitations import collect_excited_state_gradient
            collect_excited_state_gradient(self, f, indices)

        # Collect excitations and transitions
        if 'EXCITATIONS' in self.calctype:
            from .excitations import collect_excitations
            collect_excitations(self, f, indices)

        # Collect polarizability
        if 'POLARIZABILITY' in self.calctype:
            from .polarizability import collect_polarizability
            collect_polarizability(self, f, indices)

        # Collect frequencies and IR intensities
        if 'FREQUENCIES' in self.calctype:
            from .vibrations import collect_frequencies
            collect_frequencies(self, f, indices)
#        if 'DIM' in self.calctype:
#            from .dim import collect_dim
#            collect_dim(self, f, indices)

    def __collect_symmetry(self, f, indices):
        '''Collect the symmetry of the system, if it exists.'''
        if 'SYMMETRY' in indices:
            ix = indices['SYMMETRY']
            self.symmetry = f[ix].split()[-1]

    def __collect_dipole(self, f, indices):
        '''Collect the dipole moment of the system in a.u.'''
        import re
        # Find the dipole moment this way if the keyword "dipole" is 
        # inside of the property block of an NWChem input file.  This
        # is preferred if you want more digits in the dipole moment.
        from chemPackage.constants import AU2DEBYE
        dipole = re.compile(
                            r'''
                            \s*
                            (DM[X,Y,Z])
                            \s+
                            ([-0-9.]+)
                            \s+
                            DM[X,Y,Z]EFC
                            \s+
                            [-0-9.]+                          
                            ''', re.VERBOSE)
        if 'DIPOLES 1' in indices:
            s = indices['DIPOLES 1']
            e = next(i for i, x in enumerate(f[s:], s) 
                                                 if 'NWChem Input Module' in x)
            dipoles = []
            for line in f[s:e]:
                m = dipole.match(line)
                if (m):
                    dipoles.append(m.group(2))
            self.dipole = array(dipoles[0:3], dtype=float)
        # In general, the dipole moment is hidden in the output in a 
        # table near the multipole analysis.  The "L" column indicates 
        # the level of the multipole expansion (0=monopole, 1=dipole, 
        # 2=quadrupole). RHF/ROHF prints a little differently than
        # UHF for some reason, so we need two ways of matching. 
        dipole1 = re.compile(
                             r'''
                             \s*
                             ([1])
                             \s+
                             ([01])
                             \s+
                             ([01])
                             \s+
                             ([01])
                             \s+
                             ([-0-9.]+)
                             \s+
                             [-0-9.]+                          
                             \s+
                             [-0-9.]+                          
                             \s+
                             ([-0-9.]+)                         
                             ''', re.VERBOSE)
        dipole2 = re.compile(
                             r'''
                             \s*
                             ([1])
                             \s+
                             ([01])
                             \s+
                             ([01])
                             \s+
                             ([01])
                             \s+
                             ([-0-9.]+)
                             \s+
                             [-0-9.]+                          
                             \s+
                             ([-0-9.]+)                          
                             ''', re.VERBOSE)
        # The line to search for.
        if 'DIPOLES 2' in indices:
            s = indices['DIPOLES 2']
            e = s + 11
            dipoles = []
            dipoles_nuc = []
            for line in f[s:e]:
                m = dipole1.match(line)
                m1 = dipole2.match(line)
                if (m):
                    # For DFT and UHF
                    dipoles.append(m.group(5))
                    dipoles_nuc.append(m.group(6))
                elif (m1):
                    # For RHF and ROHF
                    dipoles.append(m1.group(5))
                    dipoles.append(m1.group(6))
            self.dipole = array(dipoles[0:3], dtype=float)
            self.dipoles_nuc = array(dipoles_nuc[0:3], dtype=float)
        self.dipole = AU2DEBYE(self.dipole)
#        self.dipoles_nuc = nuc*AU2DEBYE(self.dipoles)
        self.dipoles_nuc = AU2DEBYE(self.dipoles_nuc)

    def __collect_charges(self, f, indices):
        '''Collect the atomic charges'''
        if 'CHARGES' in indices:
            # Collect the Mulliken population analysis.  It is put in terms of
            # the total number of electrons - we convert it to be the number of
            # electrons added or removed compared to the free atom.
            s = indices['CHARGES']
            e = next(i for i, x in enumerate(f[s:], s) if not x.strip())
            atmnum = array([x.split()[2] for x in f[s:e]], dtype=float)
            chrg   = array([x.split()[3] for x in f[s:e]], dtype=float)
            self.charges.mulliken = atmnum - chrg
        
    def __collect_energy(self, f, indices):
        '''Collect the energy of the system'''

        # Reference table for conversion of many-body perturbation 
        # theory types to simpler names.
        mbpt_convert = {'MBPT(2)' : 'MP2', 'MBPT(3)' : 'MP3', 
                        'MBPT(4)' : 'MP4',}

        # The ground state reference can be HF or DFT.

        if 'DFT' in self.calctype:
            self.energy = {}
            if 'RESTRICTED' in self.calctype:
                a = indices['RESTRICTED DFT']
            elif 'UNRESTRICTED' in self.calctype:
                a = indices['UNRESTRICTED DFT ALPHA']
            ix = next(i for i, x in renumerate(f[:a], a)
                                        if 'Total iterative time' in x)
            ix = ix - 9
            self.energy['total DFT'] = float(f[ix].split('=')[-1].strip())
            self.energy['one-electron'] = float(f[ix+1].split('=')[-1].strip())
            self.energy['Coulomb'] = float(f[ix+2].split('=')[-1].strip())
            self.energy['XC'] = float(f[ix+3].split('=')[-1].strip())
            self.energy['nuc. repulsion'] = float(f[ix+4].split('=')[-1].strip())       
            self.energy['total'] = self.energy['total DFT']
        elif 'HF' in self.calctype:
            if 'RESTRICTED' in self.calctype:
                ix = indices['RESTRICTED HF ENERGIES']
            elif 'UNRESTRICTED' in self.calctype:
                ix = indices['UNRESTRICTED HF ENERGIES']
            self.energy = {}
            self.energy['total HF'] = float(f[ix].split('=')[-1].strip())
            self.energy['one-electron'] = float(f[ix+1].split('=')[-1].strip())
            self.energy['two-electron'] = float(f[ix+2].split('=')[-1].strip())
            self.energy['nuc. repulsion'] = float(f[ix+3].split('=')[-1].strip())
        if 'DIM' in self.calctype:
            #ix = indices['DIM/QM ENERGY'] + 2
            ix = indices['DIM/QM ENERGY'] + 1
            #self.energy['DIM el']  = float(f[ix].split('=')[-1].strip())
            #self.energy['DIM nuc'] = float(f[ix+1].split('=')[-1].strip())
            #self.energy['DIM total'] = float(f[ix+2].split('=')[-1].strip())
            self.energy['DIM total'] = float(f[ix].split(' ')[-1].strip())
            #self.dimqm_energy = self.energy['DIM total']
            self.dimqm_energy = self.energy['DIM total']

        # For MP2 and coupled cluster calculations submitted through
        # keyword blocks and not TCE.  These are mostly for energy
        # calculations since TCE is the most versatile method for
        # building highly correlated wavefunctions in NWChem.

        if 'MP2' in self.calctype and not 'TCE' in self.calctype: 
            # Index to start the search
            s = indices['MP2 BLOCK ENERGY']
            # Find where the MP2 energy is
            arr = [i for i,x in enumerate(f[s:],s) if 'Correlation energy' in x]
            # Collect MP2 energy information
            for ix in arr:
                self.energy['MP2 correlation'] = float(f[ix].split()[2])
                self.energy['MP2 singlet'] = float(f[ix+1].split()[2])
                self.energy['MP2 triplet'] = float(f[ix+2].split()[2])
                self.energy['total MP2'] = float(f[ix+3].split()[3])
        elif 'CCSD' in self.calctype and not 'TCE' in self.calctype:
            # MP2 correlation and total energy.
            ix = indices['CCSD BLOCK MP2 ENERGY']
            self.energy['MP2 correlation'] = float(f[ix].split()[3])
            self.energy['total MP2'] = float(f[ix+1].split()[3])
            # CCSD correlation and total energy.
            ix = indices['CCSD BLOCK ENERGY']
            self.energy['CCSD correlation'] = float(f[ix].split()[3])
            self.energy['total CCSD'] = float(f[ix+1].split()[3])

        # For TCE, the higher-level correlation method can be:  
        # configuration interaction, many-body perturbation theory, 
        # and coupled cluster.  

        if 'TCE' in self.calctype:
            s = indices['TCE HEAD']

            for ln in f[s:]:
                if 'correlation' in ln:
                    if ln.split()[0] in mbpt_convert.keys():
                        theory = mbpt_convert[ln.split()[0]]
                    else:
                        theory = ln.split()[0]
                    theorystr = theory + ' correlation'
                    self.energy[theorystr] = float(ln.split()[6])
                if 'total' in ln and 'bytes' not in ln:
                    if 'total' == ln.split()[1]:
                        if ln.split()[0] in mbpt_convert.keys():
                            theory = mbpt_convert[ln.split()[0]]
                        else:
                            theory = ln.split()[0]
                        theorystr = 'total ' + theory
                        self.energy[theorystr] = float(ln.split()[6])

    def __det_calc_type(self):
        '''The calculation type of this file is determined based on the
        entries found in key and subkey.  This is done for cross-package
        support.
        '''
        if 'TASK' in self.key:
            for k in self.key['TASK']:
                if 'MCSCF' in k.upper():
                    self.calctype.add('MCSCF')
                elif 'SCF' in k.upper():
                    self.calctype.add('HF')
                if 'DFT' in k.upper(): self.calctype.add('DFT')
                if 'CCSD' in k.upper(): self.calctype.add('CCSD') # For CCSD block
                if 'MP2' in k.upper(): self.calctype.add('MP2')   # For MP2 block
                if 'TCE' in k.upper():
                    self.calctype.add('TCE')
                    if 'CCSD' in self.subkey:
                        self.calctype.add('CCSD')
                    elif ('MBPT2' in self.subkey or
                          'MP2' in self.subkey):
                        self.calctype.add('MP2')
                    elif 'CISD' in self.subkey:
                        self.calctype.add('CISD')
                    # Calculations may be done with a HF or DFT 
                    # ground state reference.
                    if 'SCF' in self.key:
                        self.calctype.add('HF')
                    else:
                        self.calctype.add('DFT')
                if 'TDDFT ENERGY' in k.upper():
                    self.calctype.add('EXCITATIONS')
                    if 'CDSPECTRUM' in self.subkey:
                        self.calctype.add('CD SPECTRUM')
                    if 'VELOCITY' in self.subkey:
                        self.calctype.add('VELOCITY')
                if 'OPTIMIZE' in k.upper():
                    self.calctype.add('GEOMETRY')
                    if 'TDDFT' in k.upper(): 
                        self.calctype.add('EXCITED STATE')
                        self.calctype.add('EXCITATIONS')
                if 'GRADIENT' in k.upper():
                    self.calctype.add('GEOMETRY')
                    if 'TDDFT' in k.upper(): 
                        self.calctype.add('EXCITED STATE')
                        self.calctype.add('EXCITATIONS')
                if 'FREQUENCIES' in k.upper(): self.calctype.add('FREQUENCIES')
                if 'FREQ' in k.upper(): self.calctype.add('FREQUENCIES')
        if 'AORESPONSE' in self.key or 'RESPONSE' in self.subkey:
            # NWChem calculates both polarizabilities and the
            # optical rotation Beta tensor when response/aoresponse
            # are used.
            self.calctype.add('POLARIZABILITY')
            self.calctype.add('OPTICAL ROTATION')
            if 'DAMPING' in self.subkey:
                self.calctype.add('FD')
            else:
                self.calctype.add('STATIC')
        if 'SET' in self.key:
            for k in self.key['SET']:
                if 'TCE:LINERESP T' in k.upper():
                    self.calctype.add('POLARIZABILITY')
                    if 'SET TCE:AFREQ' in k.upper():
                        self.calctype.add('FD') # Not damped
                    else:
                        self.calctype.add('STATIC')
        if 'TCE' in self.calctype:
            for k in self.key['TCE'][0]:
                if 'NROOTS' in k.upper():
                    self.calctype.add('EXCITATIONS')
        if 'ODFT' in self.subkey or 'UHF' in self.subkey:
            self.calctype.add('UNRESTRICTED')
        else:
            # Closed shell DFT, RHF, ROHF
            self.calctype.add('RESTRICTED')
        if 'DIMQM' in self.key:
            self.calctype.add('DIM')
        if 'FREQ' in self.key:
            # Need to check that isotopic substitution is not being done,
            # and that we aren't reusing the Hessian.
            for k in self.key['FREQ']:
                for item in k:
                    if 'REUSE' in item.upper():
                        self.calctype.add('REUSE HESSIAN') 
                    elif 'MASS' in item.upper():
                        self.calctype.add('ISOTOPE SUBSTITUTION')

def renumerate(sequence, start=0):
    n = start
    for elem in reversed(sequence):
        yield n, elem
        n -= 1
