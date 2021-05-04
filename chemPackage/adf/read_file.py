from __future__ import print_function, division

def read_file(self):
    '''Read in the file and store where major sections begin.'''

    # Collect all data into memory for data retention
    with open(self.filename) as fl:
        f = tuple([line.rstrip() for line in fl])

    # For an input file, grab start and end of input block and return
    if self.filetype != 'out':
        istrt = next(i for i, x in enumerate(f) if r'$ADFBIN/adf' in x)
        iend = next(i for i, x in enumerate(f[istrt:], istrt) if x == 'eor')
        indices = { 'INPUT START' : istrt, 'INPUT END' : iend }
        return f, indices
    # Otherwise, read the entire output file
    if self.project=='all':
        # Define lines that are accociated with various keys or properties.
        # Since these lines may appear more than once, we define three groups:
        # one where we want the first appearence, one where we want the last
        # appearance, and one where we want every appearence
        # The keys are the lines to find, and the first index of the value is the
        # associated propery, and the second is the number to add the the line
        # number to find the start of where to collect the property.
        #
        # Note that since we are defining the entire line, not just a part,
        # we can search these lines as part of a set which are designed for
        # searching.  If we defined a part of the line, we would have to scan each
        # line, being much slower.
        first = {
                #############
                # Input block
                #############

                # End of the input block
                'end input':
                                                                  ['INPUT END', 0],
                #############################
                # General / Common properties
                #############################

                # Total energies
                'Pauli Repulsion':
                                                                     ['ENERGY', 0],

                ##########
                # Geometry
                ##########

                # Initial geometry
                ' =====                            X Y Z'
                '                    CHARGE':
                                                           ['INITIAL GEOMETRY', 3],
                # Initial geometry for a COSMO calculation
                ' =====                            X Y Z'
                '                    CHARGE                     RADIUS':
                                                     ['INITIAL GEOMETRY COSMO', 3],

                ##############################################
                # Polarizability / Hyperpolarizability / Raman
                ##############################################

                # RESPONSE polarizability tensors
                'THE DIPOLE-DIPOLE POLARIZABILITY TENSOR:': 
                                                                   ['RESPONSE', 0],
                # Hyperpolarizability invarients
                ' Normal termination of HYPERPOL program part':
                                             ['HYPERPOLARIZABILITY INVARIENTS' ,0],
                #########################
                # Vibrational frequencies
                #########################

                # Normal modes
                ' Vibrations and Normal Modes  ***  (cartesian coordinates,'
                ' NOT mass-weighted)  ***':
                                                               ['NORMAL MODES', 5],
                # IR frequencies
                ' List of All Frequencies:':
                                                                         ['IR', 9],
                # Analytical frequencies
                '=== CALCULATE ANALYTICAL SECOND DERIVATIVES OF THE ENERGY ===':
                                                     ['ANALYTICAL FREQUENCIES', 0],

                #Dipole Derivatives
                ' Derivatives:':
                                                         ['DIPOLE DERIVATIVES', 0],

                ####################################
                # Excited State Gradients / Geometry
                ####################################

                # Excited state gradients
                '     ++Total Excitation Energy Gradient++': 
                                                    ['EXCITED STATE GRADIENTS', 5],
                #############
                # Excitations
                #############

                # End of the excitation parts
                ' Normal termination of EXCITATION program part':
                                                            ['EXCITATIONS END', 0],
                # Unrestricted excitations
                ' All SPIN-UNRESTRICTED excitation energies':
                                                          ['SPIN-UNRESTRICTED', 4],
                # Singlet-singlet excitations
                ' All SINGLET-SINGLET excitation energies': 
                                                            ['SINGLET-SINGLET', 4],
                # Singlet-triplet excitations
                ' All SINGLET-TRIPLET excitation energies': 
                                                            ['SINGLET-TRIPLET', 4],
                #####
                # DIM
                #####

                # DIM coordinates
                ' DIM-SYSTEM':
                                                            ['DIM COORDINATES', 5],
                # DIM dipole moment
                ' D I M / Q M   R E S U L T S':
                                                          ['DIM DIPOLE MOMENT', 0],
                # DIM/QM interaction energy
                ' DIM/QM Interaction Energy :':
                                                              ['DIM/QM ENERGY', 0],
                # DIM energy
                ' DIM System Energy :':
                                                                 ['DIM ENERGY', 0],
                ###########
                # Technical
                ###########

                # Exit on error
                ' ADF EXIT called': 
                                                                      ['ERROR', 1],
                # Detailed timing section
                '   User CPU   Sys CPU   Elapsed   Ncall   Timer Name':
                                                             ['TIMING VERBOSE', 2],
                # Parallelization
                ' Rank   Node Name                           '
                '   NodeID   MyNodeRank  NodeMaster': 
                                                        ['NEW PARALLELIZATION', 0],
                # Old version of parallelization
                '  Kid Node                   Speed          Task ID   Status':
                                                        ['OLD PARALLELIZATION', 0],
                # End of timing
                ' Buffered I/O statistics':
                                                                 ['TIMING END', 0],
                }

        last =  {
                #############
                # Input block
                #############

                # Start of input file
                '(INPUT FILE)':
                                                                ['INPUT START', 0],
                #############################
                # General / Common properties
                #############################

                # Symmetry of the system
                ' S Y M M E T R Y ,   E L E C T R O N S':
                                                                   ['SYMMETRY', 3],
                # Total dipole moment
                ' Dipole Moment  ***  (Debye)  ***':
                                                                   ['DIPOLE L', 3],
                # Quadrupole Moment - Buckingham Convention
                ' Quadrupole Moment (Buckingham convention)  ***  (a.u.)  ***':
                                                        ['QUADRUPOLE MOMENT L', 4],
                # Quadrupole Moment - Fit density
                ' Represented molecular multipole moments':
                                                       ['QUADRUPOLE FIT DENS', 13],
                # Atomic multipole moments
                ' Atomic electronic multipole moments from SCF equations (a.u.)':
                                                   ['ATOMIC MULTIPOLE MOMENTS', 5],
                # Multipole-derived charges
                ' Multipole derived atomic charges (a.u.)':
                                                  ['MULTIPOLE-DERIVED CHARGES', 9],
                # Hirshfeld fragment charges
                ' H I R S H F E L D   C H A R G E   A N A L Y S I S':
                                                         ['HIRSHFELD CHARGES', 10],
                # Mulliken charges
                ' M U L L I K E N   P O P U L A T I O N S':
                                                          ['MULLIKEN CHARGES', 10],
                # Voronoi chages
                '  Atom                Initial                    OrthFrag                   SCF':
                                                            ['VORONOI CHARGES', 3],
                # SFO (Symmetrized Fragment Orbitals) Basis
                '  indx  incl.CFs)   Occup   Orb.Energy   FragmentType  Coeff.   '
                'Orbital     on Fragment':
                                                                  ['SFO BASIS', 2],

                # Cartesian STO basis
                ' (power of) X  Y  Z  R     Alpha  on Atom':
                                                                  ['STO BASIS', 3],
                # FDE coupled polarizability
                '                  X             Y             Z':
                                                        ['FDEC POLARIZABILITY', 2],

                ####################################
                # Excited State Gradients / Geometry
                ####################################

                # Excited state dipole moment
                ' *** Finished Excitation Energy Gradients ***': 
                                                      ['EXCITED STATE DIPOLE', -2],
                ###########
                # Technical
                ###########

                # Normal end
                '                             A D F   E X I T':
                                                                   ['ADF EXIT', 1],
                }

        each =  {
                #############
                # Input block
                #############

                # Alternative input block ending if the first doesn't work
                # First of these after the input start
                ' ******************************************'
                '*************************************':
                                                          ['INPUT END BACKUP', -2],
                #############################
                # General / Common properties
                #############################

                # Total dipole moment
                ' Dipole Moment  ***  (Debye)  ***': 
                                                                     ['DIPOLE', 3],

                # Quadrupole Moment - Buckingham Convention
                ' Quadrupole Moment (Buckingham convention)  ***  (a.u.)  ***':
                                                          ['QUADRUPOLE MOMENT', 4],

                # Orbital block
                '       E(eV)  Occ       MO           %     SFO (first member)'
                '   E(eV)  Occ   Fragment':
                                                                   ['ORBITALS', 3],
                ##########
                # Geometry
                ##########

                # Optimized geometry
                '                   X           Y           Z              '
                'X           Y           Z       (0:frozen, *:LT par.)':
                                                         ['OPTIMIZED GEOMETRY', 2],
                ##############################################
                # Polarizability / Hyperpolarizability / Raman
                ##############################################

                # AORESPONSE polarizability
                ' POLARIZABILITY':
                                                                 ['AORESPONSE', 0],
                # AORESPONSE Optical Rotation
                ' OPTICAL ROTATION':
                                                           ['OPTICAL ROTATION', 0],
                # Polarizability derivatives in a Raman calculation
                '  Polarizability Derivatives':
                                                 ['POLARIZABILITY DERIVATIVES', 0],
                # Hyperpolarizability tensors
                '     The STATIC hyperpolarizability tensor beta':
                                                         ['HYPERPOLARIZABILITY',0],
                # Hyperpolarizability aoresponse tensors
                ' beta                    real         imaginary':
                                                         ['HYPERPOLARIZABILITY',1],
                # Second hyperpolarizability aoresponse tensors
                ' gamma                             real              imaginary':
                                                  ['SECOND HYPERPOLARIZABILITY',1],
                # Frequency-dependent Raman scattering factors
                '  Raman Scattering Factors':
                                                                   ['FD RAMAN', 9],
                '                Frequency          Depolarization ratio'
                # Static Raman scattering factors
                '     Raman Intensity (deg. not counted)':
                                                               ['STATIC RAMAN', 3],
                # VROA intensities
                ' VROA Intensities  (10E3 Ang**4/amu)':
                                                           ['VROA INTENSITIES', 9],
                # Dipole-Quadrupole Polarizabilities
                ' DIPOLE-QUADRUPOLE POLARIZABILITY':
                                                                   ['A-TENSOR', 8],
                # Polarizability from an AORESPONSE calculation
                '      DIPOLE-QUADRUPOLE POLARIZABILITY':
                                                                   ['A-TENSOR', 2],
                # Dipole-Octupole polarizability
                ' DIPOLE-OCTUPOLE POLARIZABILITY':
                                                                   ['B-TENSOR', 8],
                # HYPERPOLARIZABILITY dipole-dipole-dipole
                ' HYPERPOLARIZABILITY dipole-dipole-dipole':
                                                    ['HYPERPOLARIZABILITY_DDD', 4],
                # HYPERPOLARIZABILITY dipole-dipole-quadrupole
                ' HYPERPOLARIZABILITY dipole-dipole-quadrupole':
                                                    ['HYPERPOLARIZABILITY_DDQ', 4],
                # HYPERPOLARIZABILITY dipole-quadrupole-dipole
                ' HYPERPOLARIZABILITY dipole-quadrupole-dipole':
                                                    ['HYPERPOLARIZABILITY_DQD', 4],
                # HYPERPOLARIZABILITY quadrupole-dipole-dipole
                ' HYPERPOLARIZABILITY quadrupole-dipole-dipole':
                                                    ['HYPERPOLARIZABILITY_QDD', 4],
                # HYPERPOLARIZABILITY dipole-quadrupole-quadrupole
                ' HYPERPOLARIZABILITY dipole-quadrupole-quadrupole':
                                                    ['HYPERPOLARIZABILITY_DQQ', 4],
                # HYPERPOLARIZABILITY quadrupole-dipole-quadrupole
                ' HYPERPOLARIZABILITY quadrupole-dipole-quadrupole':
                                                    ['HYPERPOLARIZABILITY_QDQ', 4],
                # HYPERPOLARIZABILITY quadrupole-quadrupole-dipole
                ' HYPERPOLARIZABILITY quadrupole-quadrupole-dipole':
                                                    ['HYPERPOLARIZABILITY_QQD', 4],
                # HYPERPOLARIZABILITY quadrupole-quadrupole-quadrupole
                ' HYPERPOLARIZABILITY quadrupole-quadrupole-quadrupole':
                                                    ['HYPERPOLARIZABILITY_QQQ', 4],
                # Polarizability derivatives more easily printed
                '  Polarizability (a) Derivatives':
                                               ['POLARIZABILITY DERIVATIVES R', 3],
                '  Imaginary Polarizability (a) Derivatives':
                                               ['POLARIZABILITY DERIVATIVES I', 3],
                # Optical rotation derivatives
                '  Optical activity (G) Derivatives':
                                             ['OPTICAL ROTATION DERIVATIVES R', 3],
                '  Imaginary Optical activity (G) Derivatives':
                                             ['OPTICAL ROTATION DERIVATIVES I', 3],
                # A-tensor derivatives
                '  Dipole-Quadrupole Polarizability (A) Derivatives':
                                                     ['A-TENSOR DERIVATIVES R', 3],
                '  Imaginary Dipole-Quadrupole Polarizability (A) Derivatives':
                                                     ['A-TENSOR DERIVATIVES I', 3],
                # Quad-Quad polarizability using new AORESPONSE routine
                ' QUADRUPOLE-QUADRUPOLE POLARIZABILITY':
                                                                   ['C-TENSOR', 8],
                # Quad-Mag dipole polarizability (Dtensor)
                ' QUADRUPOLE-MAGNETIC DIPOLE Polarizability tensor:':
                                                                   ['D-TENSOR', 2],
                # Magnetic dipole-magnetic dipole polarizability (magnetizability)
                ' MAGNETIZABILITY':
                                                            ['MAGNETIZABILITY', 6],
                #############
                # Excitations
                #############

                # Transition dipole moments
                ' no.  E/eV          f                       mu (x,y,z)':
                                                                        ['TDM', 2],
                # Excitations energies calculated with exact algorithm
                '    Nr          Excitation energy             '
                'Oscillator Strength':
                                                          ['EXACT EXCITATIONS', 2],
                # Excitations energies calculated with Davidson algorithm
                ' no.  E/a.u.        E/eV      f           dE/a.u.':
                                                       ['DAVIDSON EXCITATIONS', 2],
                # Transitions in the excitations
                ' Major MO -> MO transitions for the above excitations':
                                                                ['TRANSITIONS', 5],
                #####
                # DIM
                #####

                # DIM dipoles and charges
                ' Induced dipoles and charges for each DIM atom :':
                                                    ['DIM DIPOLES AND CHARGES', 0],
                # DIM dipoles
                ' Induced dipoles for each DIM atom :':
                                                                ['DIM DIPOLES', 0],
                # Old DIM dipoles and charges
                '   Induced dipoles and charges for each DIM/QM atom':
                                                                ['OLD DIM D&C', 0],
                # DIM FD polarizability, real
                ' Polarizability tensor for DIM system: Frequency-Dependent Real':
                                                 ['DIM FD REAL POLARIZABILITY', 0],
                # DIM FD polarizability, imag
                ' Polarizability tensor for DIM system: Frequency-Dependent Imag':
                                                 ['DIM FD IMAG POLARIZABILITY', 0],
                # DIM static polarizability
                ' Polarizability tensor for DIM system: Static':
                                                  ['DIM STATIC POLARIZABILITY', 0],

        
                # DIM Total Dipole of the system. real
                ' Real Part':
                                                 ['DIM DIPOLE REAL', 1],
     
                ' Imaginary Part':
                                                 ['DIM DIPOLE IMAG', 1],


                ##########################
                # Frozen Density Embedding
                ##########################
                # Embedding Energies
                'Total electrostatic':             ['FDE ELECTROSTATIC ENERGY', 0],
                '                |Error in J|':                   ['FDE ERROR', 0],
                'Total interaction energy':          ['FDE INTERACTION ENERGY', 0],
                'Subsystem A total energy':                 ['FDE SUBSYSTEM A', 0],
                'Subsystem B total energy':                 ['FDE SUBSYSTEM B', 0],
                'Complete subsystem DFT energy':         ['FDE ALL SUBSYSTEMS', 0],
                'Complete embedding energy':           ['FDE EMBEDDING ENERGY', 0],
                # Kinetic energies
                'KSCED integrals: T[1]                    =':
                                                                    ['FDE T A', 0],
                'KSCED integrals: T[rho2]          (appr) :':
                                                                    ['FDE T B', 0],
                'KSCED integrals: T[rho1+rho2]     (appr) :':
                                                                   ['FDE T AB', 0],
                'KSCED ENERGY: T[1+2]-T[2]                =':
                                                               ['FDE T AB - B', 0],
                'KSCED ENERGY: Tnad[1,2]                  =':
                                                           ['FDE T AB - A - B', 0],
                # Kinetic energies (LDA)
                'KSCED integrals: T[1](lda)               =':
                                                                ['FDE T A LDA', 0],
                'KSCED integrals: T[rho2](lda)     (appr) :':
                                                                ['FDE T B LDA', 0],
                'KSCED integrals: T[rho1+rho2](lda)(appr) :':
                                                               ['FDE T AB LDA', 0],
                'KSCED ENERGY: T[1+2]-T[2](lda)           =':
                                                           ['FDE T AB - B LDA', 0],
                'KSCED integrals: Tnad[1,2](lda)          =':
                                                       ['FDE T AB - A - B LDA', 0],
                # Potentials
                'KSCED integrals: int(rho1*vtnad)         =':
                                                                     ['FDE VT', 0],
                'KSCED integrals: int(rho1*vtnad(lda))    =':
                                                                 ['FDE VT LDA', 0],
                'KSCED ENERGY: Enuc2[1]                   =':
                                                                   ['FDE ENUC', 0],
                'KSCED ENERGY: J[1,2] (2=fitted)          =':
                                                                 ['FDE J[A,B]', 0],
                'KSCED ENERGY: Enuc2[1]+J[1,2] (sum)      =':
                                                          ['FDE ENUC + J[A,B]', 0],
                # XC energies
                'KSCED integrals: Exc[rho1]               =':
                                                                   ['FDE XC A', 0],
                'KSCED integrals: Exc[rho2]        (appr) :':
                                                                   ['FDE XC B', 0],
                'KSCED integrals: Exc[rho1+rho2]   (appr) :':
                                                                  ['FDE XC AB', 0],
                'KSCED ENERGY: Exc[1+2]-Exc[2]            =':
                                                              ['FDE XC AB - B', 0],
                'KSCED ENERGY: Excnad[1,2]                =':
                                                          ['FDE XC AB - A - B', 0],
                # XC energies (LDA)
                'KSCED integrals: Exc[rho1](lda)          =':
                                                               ['FDE XC A LDA', 0],
                'KSCED integrals: Exc[rho2](lda)   (appr) :':
                                                               ['FDE XC B LDA', 0],
                'KSCED integrals: Exc[rho1+rho2]lda(appr) :':
                                                              ['FDE XC AB LDA', 0],
                'KSCED ENERGY: Exc[1+2]-Exc[2](lda)       =':
                                                          ['FDE XC AB - B LDA', 0],
                'KSCED ENERGY: Exc[1+2]-Exc[1]-Exc[2](lda)=':
                                                      ['FDE XC AB - A - B LDA', 0],
                # coupled excitation energies
                '                 List of coupled states':
                                                                ['FDEU STATES', 5],
                '  Coupled Subsystem Excitation energies:':
                                                              ['FDEC ENERGIES', 4],
                '  Electric transition dipole moments:':
                                                        ['FDEC DIPOLE MOMENTS', 4],

                ###########
                # Technical
                ###########

                # Timing section.  First after ADF EXIT
                ' Timing Statistics': 
                                                                     ['TIMING', 0],
                }

    elif self.project=='SERS':
        # Only account for what's needed for SERS and nothing else. This saves lots of time
        # Such method of re-grouping keywords may be applied to other projects, e.g. hpol.
        # -- Pengchong, Nov. 2016
        # Zhongwei: include the 'SEHRS' part within the 'SERS' project
        first = {
                #############
                # Input block
                #############

                # End of the input block
                'end input':
                                                                  ['INPUT END', 0],
                #############################
                # General / Common properties
                #############################

                # Total energies
                'Pauli Repulsion':
                                                                     ['ENERGY', 0],

                ##########
                # Geometry
                ##########

                # Initial geometry
                ' =====                            X Y Z'
                '                    CHARGE':
                                                           ['INITIAL GEOMETRY', 3],
                # Initial geometry for a COSMO calculation
                ' =====                            X Y Z'
                '                    CHARGE                     RADIUS':
                                                     ['INITIAL GEOMETRY COSMO', 3],

                #########################
                # Vibrational frequencies
                #########################

                # Normal modes
                ' Vibrations and Normal Modes  ***  (cartesian coordinates,'
                ' NOT mass-weighted)  ***':
                                                               ['NORMAL MODES', 5],
                # Analytical frequencies
                '=== CALCULATE ANALYTICAL SECOND DERIVATIVES OF THE ENERGY ===':
                                                     ['ANALYTICAL FREQUENCIES', 0],

                #Dipole Derivatives
                ' Derivatives:':
                                                         ['DIPOLE DERIVATIVES', 0],

                }

        last =  {
                #############
                # Input block
                #############

                # Start of input file
                '(INPUT FILE)':
                                                                ['INPUT START', 0],
                }
        each =  {
                ##############################################
                # Polarizability / Hyperpolarizability / Raman
                ##############################################

                # AORESPONSE polarizability
                ' POLARIZABILITY':
                                                                 ['AORESPONSE', 0],
                # Hyperpolarizability aoresponse tensors
                ' beta                    real         imaginary':
                                                         ['HYPERPOLARIZABILITY',1],
                # Second hyperpolarizability aoresponse tensors
                ' gamma                             real              imaginary':
                                                  ['SECOND HYPERPOLARIZABILITY',1],
                }
         
    # Loop over the file and find the index where each line appears
    #search_lines = set(first.keys()+last.keys()+each.keys())
    search_lines = set(sorted(first)+sorted(last)+sorted(each))
    indices = {}
    for i, line in enumerate(f):
        # Use the broad search in a set for fastest searching.
        if ' SCF MODERATELY CONVERGED' == line:
            print(self.filename, 'WARNING: SCF MODERATELY CONVERGED')
        if line in search_lines:
            # If the line is in the first dict, store then remove from search
            if line in first:
                indices[first[line][0]] = i+first[line][1]
                del first[line]
                search_lines.remove(line)
            # If the line is in a last dict, store and overwrite
            elif line in last:
                indices[last[line][0]] = i+last[line][1]
            # Otherwise, append this value to the index
            else:
                try:
                    indices[each[line][0]].append(i+each[line][1])
                except KeyError:
                    indices[each[line][0]] = []
                    indices[each[line][0]].append(i+each[line][1])
        else:
            # FDE lines involve partial line matches
            for word in each.keys():
                if line[:len(word)] == word and each[word][0][:3] == 'FDE':
                    try:
                        indices[each[word][0]].append(i+each[word][1])
                    except KeyError:
                        indices[each[word][0]] = []
                        indices[each[word][0]].append(i+each[word][1])

    # Make sure that we take the timing stats immeadiately after the EXIT
    if 'ADF EXIT' in indices:
        it = indices['TIMING']
        indices['TIMING'] = next(x for x in it if x > indices['ADF EXIT'])
    else:
        if 'TIMING' in indices:
            del indices['TIMING']

    # In ADF 2014, the regular "INPUT START" keyword is not printed
    # determine the start of the input block from the "INPUT END" index
    if not 'INPUT START' in indices and 'INPUT END' in indices:
        for x in range(indices['INPUT END'],0,-1):
            if (f[x]==' =============================================================================='
                or f[x]==' Parallel Execution: Process Information'
                or '>  END' in f[x]):
                indices['INPUT START'] = x
                break

    # If no input end, take the backup after the input start
    if 'INPUT START' in indices and 'INPUT END' not in indices:
        ie, ib, ins = 'INPUT END', 'INPUT END BACKUP', 'INPUT START'
        indices[ie] = next(x for x in indices[ib] if x > indices[ins])

    #Overwrite the DIPOLE keyword for frequencies calculations
    if 'DIPOLE DERIVATIVES' in indices:
        try:
            for x in range(len(indices['DIPOLE'])-2):
                del indices['DIPOLE'][x]
            indices['DIPOLE'] = indices['DIPOLE DERIVATIVES'] -1
        except KeyError:
            pass
    else:
        # grab only the last printed dipole and quadrupole moments
        if 'DIPOLE L' in indices:
            indices['DIPOLE'] = indices['DIPOLE L']
        if 'QUADRUPOLE MOMENT L' in indices:
            indices['QUADRUPOLE MOMENT'] = indices['QUADRUPOLE MOMENT L']

    return f, indices
