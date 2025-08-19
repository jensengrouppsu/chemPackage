from ..chemdata import ChemData
from ..errorclass import CollectionError
from numpy import array, zeros, transpose


# Create the AMS class to grab and hold all the information in the file
class AMS(ChemData):
    '''Read and contain data from AMS files
    
    This class is used to read all information from either an AMS input or
    output file.  It stores this data in formats suitable for processing.  

    The only required argument for collection is the name of the file.
 
    The data available is listed below.  If the data did not exist in the
    AMS file, then the variable will default to None (or an empty set for
    subkey and calctype, or empty dict for key).
    
    '''
    # Instantiate the AMS class
    def __init__(self, name):
        '''Set up the AMS class.  The only argument is the filename.'''
        import os

        ChemData.__init__(self)
        self.program = 'AMS'
        self.project='all'
    
        # Find extention
        ftype = os.path.splitext(name)[1]
        if ftype not in ('.out', '.run', '.inp'):
            raise ValueError (ftype+' not a recognized extention')
        self.filetype = ftype[1:]
        self.filename = name

        # These keys may or may not be good
        # Create tuples of general keywords
        # TODO: Enginekeys to be collected separately, to be implemented
        # self.enginekeys = ('GEOMETRY')
        self.mixedkeys = ('ATOMS', 'EFIELD', 'FRAGMENTS', 'INTEGRATION')
        self.blockkeys = ('ANALYTICALFREQ', 'AORESPONSE', 'BASIS',
                          'CONSTRAINTS', 'DIMQM', 'DIMPAR', 'EXCITATIONS',
                          'EXCITEDGO', 'EXTERNALS', 'FDE', 'GEOVAR',
                          'GUIBONDS', 'MODIFYEXCITATION','RESPONSE', 'REMOVEFRAGORBITALS', 'SCF',
                          'SOLVATION', 'UNITS', 'VIBRON', 'XC', 'SUBEXCI', 'SUBRESPONSE')
        self.linekeys = ('A1FIT', 'BONDORDER', 'CHARGE', 'CREATE',
                         'NOPRINT', 'PRINT', 'RELATIVISTIC',
                         'SAVE', 'SCANFREQ', 'SYMMETRY', 'THERMO', 'TITLE',
                         'DEPENDENCY','DFTB RESOURCES DIR', 'BACKEND')
        self.singlekeys = ('ALLPOINTS', 'BADER', 'FORCEALDA', 'NEWDIIS',
                           'UNRESTRICTED', 'DIFFUSE', 'EXACTDENSITY',
                           'IGNOREOVERLAP', 'AOMAT2FILE', 'STOFIT', 'TOTALENERGY',
                           'AOMAT2FILE', 'ALLOWPARTIALSUPERFRAGS')

    def _collect(self, abort=False):
        '''Collect the AMS data from file.  'abort' will cause collection
        to stop on an error, instead of ignoring the error. '''

        # Save the abort status
        self._abort = abort

        # Conventions used in _collect and helper functions:
        # fi = input block
        # fo = output information block
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
        if self.project=='all':
            # Read in file
            from .read_file import read_file
            f, indices = read_file(self)
    
            # Collect input block
            for prefix in ("ADF", "DFTB", "ML"):
                start_key= f"{prefix} START"
                if start_key in indices:
                    from .input_block import collect_input
                    collect_input(self, f, indices)
            
            # Determine calculation type
            self.__det_calc_type()

            if self.filetype != 'out': return
            
            if 'DIM' in self.calctype:
                from .dim import collect_dim
                collect_dim(self, f, indices)

            if 'INITIAL GEOMETRY' in indices:
                from .coordinates import collect_geometry
                collect_geometry(self, f, indices)
    
            if 'OPTIMIZED GEOMETRY' in indices:
                from .coordinates import collect_optimized_geometry
                collect_optimized_geometry(self, f, indices)

            # Collect frequencies and IR intensities
            if "NORMAL MODES" in indices:
                from .vibrations import collect_frequencies
                collect_frequencies(self, f, indices)
                self.calctype.add('FREQUENCIES')

            if "MBH" in indices:
                from .vibrations import collect_frequencies_mbh
                collect_frequencies_mbh(self, f, indices)
                self.calctype.add('FREQUENCIES')

            # Collect polarizability
            if 'AORESPONSE' in indices:
                from .polarizability import collect_polarizability
                collect_polarizability(self, f, indices)

            # Collect G-Tensor
            if 'OPTICAL ROTATION' in indices:
                from .polarizability import collect_opticalrotation
                collect_opticalrotation(self, f, indices)
            
            # Collect A-Tensor
            if 'A-TENSOR' in indices:
                from .polarizability import collect_atensor
                collect_atensor(self, f, indices)

            # Collect C-Tensor 
            if 'C-TENSOR' in indices:
                from .polarizability import collect_ctensor
                collect_ctensor(self, f, indices)

            # Collect hyperpolarizability
            if 'HYPERPOLARIZABILITY' in indices:
                if 'RESPONSE' in self.key:
                    pass
                    # from .polarizability import collect_hyperpolarizability
                    # collect_hyperpolarizability(self, f, indices)
                else:
                    from .polarizability import collect_aoresponse_hyperpolarizability
                    collect_aoresponse_hyperpolarizability(self, f, indices)

                    # from .polarizability import collect_quadrupole_beta
                    # collect_quadrupole_beta(self, f, indices)

            # Collecct magnetizability
            if 'MAGNETIZABILITY' in indices:
                from .polarizability import collect_magnetizability
                collect_magnetizability(self, f, indices)

            self.__collect_energy(f, indices)
            # collect linear response
            if "LINEAR RESPONSE" in indices:
                from .polarizability import collect_linearresponse
                collect_linearresponse(self, f, indices)

            # Collect excitation and transtions
            if "EXCITATIONS" in self.calctype:
                if 'EXCITED STATE' in self.calctype:
                    from .excitations import collect_excited_state
                    collect_excited_state(self, f, indices)
                else:
                    from .excitations import collect_excitations
                    collect_excitations(self, f, indices)

            self.__collect_charge_analyses(f, indices)
            # Collect Hirshfeld polarizability
            self.__collect_hirshfeld_polarizability(f, indices)


    def __collect_energy(self, f, indices):
        '''Collect energy analysis.'''
        if 'ENERGY' in indices:
            ix = indices['ENERGY']
            # The energy decomposition is collected in Hartrees.  The orbial
            # analysis can be broken into symmetry groups, so collect each of
            # those.
            self.energy = {}
            fn = lambda x: float('nan') if '***' in x else float(x)
            tl = next(x for x in f[ix:] if 'Total Pauli Repulsion:' in x)
            self.energy['pauli'] = fn(tl.split()[3])
            tl = next(x for x in f[ix:] if 'Total Steric Interaction:' in x)
            self.energy['steric'] = fn(tl.split()[3])
            tl = next(x for x in f[ix:] if 'Electrostatic Interaction:' in x)
            self.energy['electrostatic'] = fn(tl.split()[2])
            tl = next(x for x in f[ix:] if 'Total Orbital Interactions:' in x)
            tp = fn(tl.split()[3])
            tl = next(x for x in f[ix:] if 'Total Bonding Energy:' in x)
            self.energy['total'] = fn(tl.split()[3])
            s = next(i for i, x in enumerate(f[ix:], ix)
                                            if x == 'Orbital Interactions') + 1
            e = next(i for i, x in enumerate(f[s:], s) if '------------' in x)
            # We will create a dict where the key is the sym group (no colon)
            # Split on the colon
            ar = [x.split(':') for x in f[s:e]]
            # First part is key.  First number in second part is value
            self.energy['orbital'] = dict(
                                [(x.strip(), fn(y.split()[0])) for x, y in ar])
            # Change HF component
            if '(Hybrid part) HF exchange' in self.energy['orbital']:
                self.energy['orbital']['HF exchange'] = (
                           self.energy['orbital']['(Hybrid part) HF exchange'])
                del self.energy['orbital']['(Hybrid part) HF exchange']
            self.energy['orbital']['total'] = tp


    def __det_calc_type(self):
        if not self.key:
            from os.path import splitext
            from chemPackage import collect, CollectionError
            from textwrap import dedent
            for ext in ('.run', '.inp'):
                inp = splitext(self.filename)[0]+ext
                try:
                    inp = collect(inp)
                except IOError:
                    pass
                else:
                    break
            else:
                from textwrap import dedent
                raise CollectionError (dedent('''\
                The chem package requires the output file to have a copy of
                the input block or have the input file be in the same directory
                as the output file. The input block is used to determine the
                calculation type and check for errors.
                '''))
            # Copy the input file's input block into this one
            self.key    = inp.key
            self.subkey = inp.subkey

        # Determine calculation type
        if 'EXCITEDGO' in self.key:
            self.calctype.add('EXCITED STATE')
        if 'EXCITATIONS' in self.key:
            # TODO: Since legacy ADF allows for both EXCITATION and EXCITATIONS
            # chemPackge also uses both version. Nasty workaround,
            # low-level code to be corrected.
            self.calctype.add('EXCITATIONS')
            self.calctype.add('EXCITATION')
        if 'RAMAN' in self.subkey:
            self.calctype.add('RAMAN')
        if 'VROA' in self.subkey:
            self.calctype.add('VROA')
            self.calctype.add('OPTICAL ROTATION')
            self.calctype.add('POLARIZABILITY')
        if 'AORESPONSE' in self.key or 'RESPONSE' in self.key:
            if 'OPTICAL ROTATION' in self.subkey:
                self.calctype.add('OPTICAL ROTATION')
            else:
                self.calctype.add('POLARIZABILITY')
            if 'HYPERPOL' in self.subkey or 'TWONPLUSONE' in self.subkey or 'BETA' in self.subkey:
            #if 'HYPERPOL' in self.subkey or 'BETA' in self.subkey:
                self.calctype.add('HYPERPOLARIZABILITY')
            if 'GAMMA' in self.subkey:
                self.calctype.add('SECOND HYPERPOLARIZABILITY')
            if 'LIFETIME' in self.subkey or 'DYNAHYP' in self.subkey:
                self.calctype.add('FD')
            else:
                self.calctype.add('STATIC')
        # if 'FREQUENCIES' in self.subkey:
        #     self.calctype.add('FREQUENCIES')
        if 'GEOMETRY' in self.key:
            it = 'ITERATIONS'
            x = [x for x in self.key['GEOMETRY'][1] if it in x.upper()]
            if 'EXCITED STATE' in self.calctype:
                try:
                    i = int(x[0].split()[1])
                except IndexError:
                    self.calctype.add('GEOMETRY')
                else:
                    if i > 1: self.calctype.add('GEOMETRY')
            elif 'FREQUENCIES' not in self.subkey:
                self.calctype.add('GEOMETRY')
            elif 'ANALYTICALFREQ' in self.key:
                self.calctype.add('GEOMETRY')
        if 'DIMQM' in self.key:
            self.calctype.add('DIM')
        if 'FRAGMENTS' in self.key:
            if self.key['FRAGMENTS'][1][0].split()[1][0:3] != 't21':
                self.calctype.add('FRAGMENT ANALYSIS')
        if 'UNRESTRICTED' in self.key:
            self.calctype.add('UNRESTRICTED')
        if 'FDE' in self.key:
            self.calctype.add('FDE')
        if 'SUBEXCI' in self.key:
            self.calctype.add('FDE')
            self.calctype.add('SUBEXCI')
        if 'SUBRESPONSE' in self.key:
            self.calctype.add('FDE')
            self.calctype.add('SUBRESPONSE')
        if 'SOLVATION' in self.key:
            self.calctype.add('COSMO')

        # All ADF is DFT calculations
        self.calctype.add('DFT')

    def __collect_charge_analyses(self, f, indices):
        '''Collects the Voroni, Hirshfeld and Multipole-derived charge analyses.'''

        if 'ATOMIC MULTIPOLE MOMENTS' in indices:
            ar = indices['ATOMIC MULTIPOLE MOMENTS']
            # Zhongwei: in case of ValueError
            try:
               self.atomic_charges = array([float(f[ar+ix].split()[2]) for ix in range(self.natoms)])
               self.atomic_dipoles = array([f[ar+ix].split()[3:6] for ix in range(self.natoms)], dtype=float)
               self.atomic_quadrupoles = array([f[ar+ix].split()[6:12] for ix in range(self.natoms)], dtype=float)
            except ValueError:
               pass

        if 'HIRSHFELD CHARGES' in indices:
            ar = indices['HIRSHFELD CHARGES']
            hirshfeld_charges = zeros((self.natoms),dtype=float)
            hcget = 0
            for i in range(len(f)-ar):
                temp = f[ar+i].split()
                if len(temp) != 3: continue
                if '===' in f[ar+i]: break
                try:
                    hirshfeld_charges[hcget] = float(temp[2])
                    hcget += 1
                except ValueError:
                    continue
                except IndexError:
                    break
                if hcget >= self.natoms: break
            try:
                self.charges.hirshfeld = hirshfeld_charges[:hcget]
            except (ValueError, TypeError):
                pass

        if 'MULLIKEN CHARGES' in indices:
            ar = indices['MULLIKEN CHARGES']
            # temporary fix
            try:
                self.charges.mulliken = array([float(f[ix+ar].split()[2]) for ix in range(self.natoms)])
            except (IndexError, TypeError):
                pass

        if 'VORONOI CHARGES' in indices:
            ar = indices['VORONOI CHARGES']
            vkeys = ['Initial Sphere', 'Initial RestCell', 'Initial NetTotal',
                     'OrthFrag Sphere', 'OrthFrag RestCell', 'OrthFrag NetTotal',
                     'SCF Sphere', 'SCF RestCell', 'SCF NetTotal', 'VDD']
            self.charges.voronoi = {}
            try:
                for ik in range(len(vkeys)):
                    self.charges.voronoi.update({vkeys[ik]: array([float(f[ix+ar].split()[ik+2]) for ix in range(self.natoms)])})
            except TypeError:
                pass
    def __collect_hirshfeld_polarizability(self, f, indices):
        '''Collects the Hirshfeld local and non-local contribution to the polarizability'''
        '''Collects the Hirshfeld charges'''
        '''Array structured as [no.Atom][iDir][jDir]'''

        # Get line numbers for Hirshfeld polarizability
        ar = []
        for i in range(len(f)):
            if 'Hirshfeld fragment' in f[i]: ar.append(i+1)
        ar = array(ar)
        if len(ar) == 0: return

        # Initialize
        nfrag = int(len(ar)/3)
        lcmplx = False
        if f[ar[0]].split()[0].upper() == 'REAL':
            lcmplx = True
            ar = ar + 1
        dtype = [complex if lcmplx else float][0]
        dipoles_loc = zeros((3,nfrag,3), dtype=dtype)
        dipoles_nonloc = zeros((3,nfrag,3), dtype=dtype)
        dipoles_all = zeros((3,nfrag,3), dtype=dtype) #Xing add
        dipoles_tot = zeros((3,3), dtype=dtype)
        charges = zeros((3,nfrag), dtype=dtype)

        # Collect values
        iar = -1
        for idir in range(3):
            for ifrag in range(nfrag):
                iar += 1
                # Collect induced dipoles: local / intrinsic
                for i in range(3):
                    if lcmplx:
                        dipoles_loc[idir][ifrag][i] = ( float(f[ar[iar]+i].split()[2])
                                                + 1j*float(f[ar[iar]+i].split()[3]) )
                    else:
                        dipoles_loc[idir][ifrag][i] = float(f[ar[iar]+i].split()[2])
                # Collect induced charges
                if lcmplx:
                    charges[idir][ifrag] = ( float(f[ar[iar]+3].split()[1])
                                         + 1j*float(f[ar[iar]+3].split()[2]) )
                else:
                    charges[idir][ifrag] = float(f[ar[iar]+3].split()[1])
                # Collect induced dipoles: non-local / charge transfer
                for i in range(3): 
                    if lcmplx:
                        dipoles_nonloc[idir][ifrag][i] = ( float(f[ar[iar]+i+4].split()[2])
                                                + 1j*float(f[ar[iar]+i+4].split()[3]) )
                    else:
                        dipoles_nonloc[idir][ifrag][i] = float(f[ar[iar]+i+4].split()[2])
                #sum non-local and local up Xing
                for i in range(3):
                    if lcmplx:
                        dipoles_all[idir][ifrag][i]=(float(f[ar[iar]+i].split()[2])+float(f[ar[iar]+i+4].split()[2])
                                                    + 1j*float(f[ar[iar]+i].split()[3])+1j*float(f[ar[iar]+i+4].split()[3]) )
                    else:
                        dipoles_all[idir][ifrag][i]=float(f[ar[iar]+i].split()[2])+float(f[ar[iar]+i+4].split()[2])               

            # Collect induced dipoles: total / identical to 'Polarizability tensor'
            for i in range(3): 
                if lcmplx:
                    dipoles_tot[idir][i] = ( float(f[ar[iar]+i+8].split()[2])
                                            + 1j*float(f[ar[iar]+i+8].split()[3]) )
                else:
                    dipoles_tot[idir][i] = float(f[ar[iar]+i+8].split()[2])


        # Return induced charges and dipoles
        # Each atomic dipole is dipole[index_atom][dir_dipole][dir_perturbation] -- Pengchong Liu, Oct. 2016
        self.hirshfeld_induced_dipoles_loc = transpose(dipoles_loc,(1,2,0))
        self.hirshfeld_induced_dipoles_nonloc = transpose(dipoles_nonloc,(1,2,0))
        self.hirshfeld_induced_dipoles_tot = dipoles_tot
        self.hirshfeld_induced_charges = transpose(charges)
        self.hirsh_pol = transpose(dipoles_all,(1,2,0)) #Xing
