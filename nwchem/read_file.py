from __future__ import print_function, division

def read_file(self):
    '''Read in the file and store where major sections begin.'''

    # Collect all data into memory for data retention
    with open(self.filename) as fl:
        f = tuple([line.rstrip() for line in fl])

    # For an input file, grab start and end of input block and return
    if self.filetype != 'out':
        return f, { 'INPUT START' : 0, 'INPUT END' : len(f) }

    # Otherwise, read the entire output file

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

            # Start of the input block
            '============================== echo of input deck '
            '==============================':
                                                            ['INPUT START', 0],
            # End of the input block
            '=================================================='
            '==============================':
                                                              ['INPUT END', 1],
            #########
            # General
            #########

            '      Symmetry information':
                                                               ['SYMMETRY', 3],
            # For counting numbers of occupied/virtual alpha and beta MOs
            # HF
            '                                 NWChem SCF Module':
                                                              ['MO NUMBER', 5],
            # DFT
            '          SCF calculation type: DFT':
                                                              ['MO NUMBER', 0],
            # Head of TCE module 
            '                   NWChem Extensible Many-Electron '
            'Theory Module':
                                                               ['TCE HEAD', 0],

            ##########################
            # Total Correlation Energy
            ##########################

            # CCSD block - MP2 energy
            ' MP2 Energy (coupled cluster initial guess)':
                                                  ['CCSD BLOCK MP2 ENERGY', 3],
            ' CCSD Energy':
                                                      ['CCSD BLOCK ENERGY', 3],
            # MP2 block
            '                   NWChem MP2 Semi-direct '
            'Energy/Gradient Module':
                                                       ['MP2 BLOCK ENERGY', 0],

            ##########
            # Geometry
            ##########

            # Number of atoms
            '            XYZ format geometry':
                                                           ['ATOMS NUMBER', 2],
 
            #############
            # Excitations
            #############

            # Head of DFT excitations
            '                                NWChem TDDFT Module':
                                                        ['DFT EXCITATIONS', 0],
            # TCE:CCSD Excitations
            ' EOM-CCSD transition moments / hartree':
                                                       ['CCSD EXCITATIONS', 0],
            
            ########
            # DIM/QM
            ########
            
            # DIM coordinates
            ' DIM COORDINATES':
                                                        ['DIM COORDINATES', 4],            
            # DIM dipole moment
            '                                     DIM/QM Results':
                                                      ['DIM DIPOLE MOMENT', 0],
            # DIM/QM Energy
            ' DIM/QM Energy':
                                                          ['DIM/QM ENERGY', 0],
            # DIM System Energy
            ' DIM System Energy':
                                                          ['DIM System Energy', 0],

            ###########
            # Technical
            ###########

            # End of file
            '                                     CITATION':
                                                               ['FILE END', 0],
            # End of file backup
            '                                  ACKNOWLEDGEMENT':
                                                        ['FILE END BACKUP', 0],
            # Job info
            '           Job information':
                                                               ['JOB INFO', 0],
            # Memory info
            '           Memory information':
                                                               ['MEM INFO', 0],
            }

    last =  {

            ###########################
            # General/Common Properties
            ###########################

            # Charges for DFT
            '      Total Density - Mulliken Population Analysis':
                                                                ['CHARGES', 5],
            # Charges for HF
            '  Mulliken analysis of the total density':
                                                                ['CHARGES', 5],
            # Dipole 1
            '          Dipole Moment':
                                                              ['DIPOLES 1', 0],
            # Dipole 2 for DFT
            '     Multipole analysis of the density':
                                                              ['DIPOLES 2', 5],
            # Dipole 2 for RHF, ROHF and UHF
            '       Multipole analysis of the density wrt the origin':
                                                              ['DIPOLES 2', 5],
            # Contributions to the total energy in HF
            '       Final RHF  results':
                                                 ['RESTRICTED HF ENERGIES', 3],
            '       Final UHF  results':
                                               ['UNRESTRICTED HF ENERGIES', 3],

            ##########
            # Orbitals
            ##########

            # HF/SCF (ROHF = RHF for a singlet state)
            '                       ROHF Final Molecular Orbital Analysis':
                                                          ['RESTRICTED HF', 0],
            '                    UHF Final Alpha Molecular Orbital Analysis':
                                                  ['UNRESTRICTED HF ALPHA', 0],
            '                     UHF Final Beta Molecular Orbital Analysis':
                                                   ['UNRESTRICTED HF BETA', 0],
            # DFT
            '                       DFT Final Molecular Orbital Analysis':
                                                         ['RESTRICTED DFT', 0],
            '                    DFT Final Alpha Molecular Orbital Analysis':
                                                 ['UNRESTRICTED DFT ALPHA', 0],
            '                     DFT Final Beta Molecular Orbital Analysis':
                                                  ['UNRESTRICTED DFT BETA', 0],            
            
            ###################
            # Vibrations and IR
            ###################

            # IR intensities
            ' Normal Eigenvalue ||           Projected Infra Red Intensities':
                                                                     ['IR', 3],
            # Normal modes start
            '             (Projected Frequencies expressed in cm-1)':
                                                           ['NMODES START', 4],
            # Normal modes end
            ' Normal Eigenvalue ||'
            '    Projected Derivative Dipole Moments (debye/angs)':
                                                            ['NMODES END', -4],

            }

    each =  {
            ##########
            # Geometry
            ##########

            # Geometry
            '  No.       Tag          Charge          '
            'X              Y              Z':
                                                               ['GEOMETRY', 2],
            ##############################################
            # Polarizability / Hyperpolarizability / Raman
            ##############################################

            # DFT polarizability
            ' DFT Linear Response polarizability / au':
                                                                ['DFT POL', 0],
            # DFT optical rotation
            ' Magnetic Dipole Response Matrix (nonzero elements):':
                                                   ['DFT OPTICAL ROTATION', 0],
            # CCSD polarizability
            ' CCSD Linear Response polarizability / au':
                                                               ['CCSD POL', 0],
            ######
            # DIM
            ######
            
            # Total DIM induced dipole moment
            ' Total induced dipole moment in DIM system :':
                                                     ['DIM DIPOLES MOMENT', 0],
            # DIM dipoles
            ' Induced dipoles for each DIM atom :':
                                                            ['DIM DIPOLES', 0],
            # DIM FD polarizability, real
            'Polarizability tensor for DIM system: Frequency-Dependent Real':
                                             ['DIM FD REAL POLARIZABILITY', 0],
            # DIM FD polarizability, imag
            'Polarizability tensor for DIM system: Frequency-Dependent Imag':
                                             ['DIM FD IMAG POLARIZABILITY', 0],
            # DIM static polarizability
            'Polarizability tensor for DIM system: Static':
                                              ['DIM STATIC POLARIZABILITY', 0],

            ############################
            # DFT and TDDFT Optimization
            ############################

            # DFT geometry optimization start and end
            '                            NWChem DFT Gradient Module':
                                                         ['DFT GRAD START', 0],
            '                 |  Time  |  1-e(secs)   |  2-e(secs)   |':
                                                           ['DFT GRAD END', 0],
            # TDDFT geometry optimization start and end
            '                           NWChem TDDFT Gradient Module':
                                                       ['TDDFT GRAD START', 0],
            ' done tddft_grad_finalize':
                                                         ['TDDFT GRAD END', 0],

            # TDDFT numerical gradients
            '                         TDDFT ENERGY GRADIENTS':
                                                       ['TDDFT GRAD START', 0],

            }

    # Loop over the file and find the index where each line appears
    search_lines = set(sorted(first)+sorted(last)+sorted(each))
    indices = {}
    for i, line in enumerate(f):
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

    # If the file end backup is used, replace file end
    if 'FILE END BACKUP' in indices:
        indices['FILE END'] = indices['FILE END BACKUP']
        del indices['FILE END BACKUP']
    # If neither were given, then call it the end of file
    elif 'FILE END' not in indices and 'FILE END BACKUP' not in indices:
        indices['FILE END'] = len(f) - 1

    # The total time is printed after the citations
    try:
        ix = indices['FILE END']
        en = enumerate
        indices['TIMING'] = next(i for i, x in en(f[ix:], ix) if 'Total' in x)
    except StopIteration:
        pass

    return f, indices
