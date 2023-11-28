
def read_file(self):
    '''Read in the file and store where major sections begin.'''

    # Collect all data into memory for data retention
    with open(self.filename) as fl:
        f = tuple([line.rstrip() for line in fl])

    # For an input file, grab start and end of input block and return
    if self.filetype != 'out':
        istrt = next(i for i, x in enumerate(f) if r'$AMSBIN/ams' in x)
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
                ' *** adf ***':
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
                # ' =====                            X Y Z'
                # '                    CHARGE':
                #                                            ['INITIAL GEOMETRY', 3],
            # Unified geometry initial geometry reading
                'Geometry':
                                                           ['INITIAL GEOMETRY', 5],

                #########################
                # Vibrational frequencies
                #########################

                # Normal modes
               #' Vibrations and Normal Modes  ***  (cartesian coordinates,'
               #' NOT mass-weighted)  ***':
                ' Normal Modes':
                                                               ['NORMAL MODES', 4],
                ' Statistical Thermal Analysis  ***  ideal gas assumed  ***':
                                                           ['NORMAL MODES END',-1],
                # IR frequencies
                ' Normal Mode Frequencies':
                                                                         ['IR', 4],
                # Analytical frequencies
                '=== CALCULATE ANALYTICAL SECOND DERIVATIVES OF THE ENERGY ===':
                                                     ['ANALYTICAL FREQUENCIES', 0],
                # Normal modes calculated using mobile block hessian
                ' Block Normal Modes (including rigid motions)':
                                                                        ['MBH', 2],
                '     CALCULATION RESULTS':
                                                                   ['MBH END', -4],
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

                }
        last =  {
                #############
                # Input block
                #############

                # Start of input file
                'ADF Engine Input':
                                                                ['INPUT START', 0],
                ###########
                # Technical
                ###########

                # Normal end
               #'                             A D F   E X I T':
                'AMS application finished. Exiting.':
                                                                   ['ADF EXIT', 0],
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
                ##########
                # Geometry
                ##########

                # Optimized geometry
                '  Index Symbol   x (angstrom)   y (angstrom)   z (angstrom)':
                                                         ['OPTIMIZED GEOMETRY', 1],
                ##############################################
                # Polarizability / Hyperpolarizability / Raman
                ##############################################

                # AORESPONSE polarizability
                ' POLARIZABILITY':
                                                                 ['AORESPONSE', 0],
                # AORESPONSE Optical Rotation
                ' OPTICAL ROTATION':
                                                           ['OPTICAL ROTATION', 0],
                # Dipole-Quadrupole Polarizabilities
                ' DIPOLE-QUADRUPOLE POLARIZABILITY':
                                                                   ['A-TENSOR', 8],
                # Quad-Quad polarizability using new AORESPONSE routine
                ' QUADRUPOLE-QUADRUPOLE POLARIZABILITY':
                                                                   ['C-TENSOR', 8],
                # Magnetic dipole-magnetic dipole polarizability (magnetizability)
                ' MAGNETIZABILITY':
                                                            ['MAGNETIZABILITY', 6],
                # Linear response function (
                ' CONDENSED LINEAR RESPONSE FUNCTION (MATRIX ELEMENTS)':
                                                             ['LINEAR RESPONSE',3],
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
                }

    # Loop over the file and find the index where each line appears
    search_lines = set(sorted(first)+sorted(last)+sorted(each))
    indices = {}
    for i, line in enumerate(f):
        # Use the broad search in a set for fastest searching.
        if ' SCF MODERATELY CONVERGED' == line:
            print(self.filename, 'WARNING: SCF MODERATELY CONVERGED')
        if 'Geometry optimization did NOT converge' == line:
            print(self.filename, 'WARNING: GEOMETRY NOT CONVERGED')
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
    
    return f, indices


