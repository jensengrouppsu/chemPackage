from __future__ import print_function, division
from .errorclass import ChemDataError

class Excitations_MOs(object):
    '''Extends the ChemData class with methods to manipulate excitations
    and molecular orbitals.
    '''

    def printHOMOLUMO(self, units='eV'):
        '''Pretty-print the :py:attr:`HOMO` - :py:attr:`LUMO` information 
        to standard output.
        
        '''
        from .constants import HART2EV
        import re

        conv = {'au': 1, 'eV':HART2EV}

        if not re.search(r'au|eV', units):
            raise ChemDataError ('Invalid choice for units: ' + units)

        assert (self.HOMO is not None) and (self.LUMO is not None), (
                                   'printHOMOLUMO(): No HOMO or LUMO present.')

        # Need to decide whether the HOMO and LUMO are associated with
        # alpha or beta spin electrons for unrestricted calculations, 
        # since the HOMO and LUMO attributes are tuples of tuples in 
        # that case.

        if 'UNRESTRICTED' in self.calctype:
            # First tuple is alpha spin, second tuple is beta spin
            homoalpha = self.HOMO[0][1]
            homobeta = self.HOMO[1][1]
            lumoalpha = self.LUMO[0][1]
            lumobeta = self.LUMO[1][1]
            if homoalpha > homobeta:
                ixh = 0
            elif homobeta > homoalpha:
                ixh = 1
            else:
                # If the HOMO for alpha and beta electrons are equal
                ixh = 0
            if lumoalpha < lumobeta:
                ixl = 0
            elif lumobeta < lumoalpha: 
                ixl = 1
            else:
                # If the LUMO for alpha and beta electrons are equal
                ixl = 0

        # Make the format string support the longest orbital ID.
        if 'RESTRICTED' in self.calctype:
            l = str(max([len(self.HOMO[0]), len(self.LUMO[0])]))
            fmt = '{0}: {1:<'+l+'} => {2:7.3f} {3} : {4}'
    
            print()
            print(fmt.format('HOMO', self.HOMO[0], 
                             conv[units] * self.HOMO[1], units, self.HOMO[2]))
            print(fmt.format('LUMO', self.LUMO[0], 
                             conv[units] * self.LUMO[1], units, self.LUMO[2]))
            print()
        else:
            l = str(max([len(self.HOMO[ixh][0]), len(self.LUMO[ixl][0])]))
            # Decide if the HOMO and LUMO are from alpha or beta spin electrons
            if ixh == 0:
                spinh = 'alpha'
            else:
                spinh = 'beta'
            if ixl == 0:
                spinl = 'alpha'
            else:
                spinl = 'beta'
            fmt = '{0}: {1:<'+l+'} {2} => {3:7.3f} {4} : {5}'

            print()
            print(fmt.format('HOMO', self.HOMO[ixh][0], spinh,
                             conv[units] * self.HOMO[ixh][1], units, 
                             self.HOMO[ixh][2]))
            print(fmt.format('LUMO', self.LUMO[ixl][0], spinl,
                             conv[units] * self.LUMO[ixl][1], units, 
                             self.LUMO[ixl][2]))
            print()

    def printOrbitals(self, thres=80, MOs=[], space=True, o1=1, 
                      o2=None, units='eV'):
        '''Pretty-print the :py:attr:`atomic orbitals <atomic_orbitals>` 
        contained in the MOs.

        By default, this will only print the atomic orbitals that make up
        80% of each MO, but this may be altered with the argument *thres*.

        By default, this will print out all MOs in the ChemData object.
        A truncated list of MOs may be given with the *MOs* argument.

        *space* indicates if a leading space will be printed.

        *o1* is the lowest orbital to print.
        *o2* is the highest, defaulting to the number of tensors.
        This option allows you to print a subset of orbitals collected.
        You may use these in conjunction with *MOs*.

        Raises :py:exc:`AssertionError` if no MOs are present.
        '''
        from .constants import HART2EV
        import re

        conv = {'au': 1, 'eV':HART2EV}

        if not re.search(r'au|eV', units):
            raise ChemDataError ('Invalid choice for units: ' + units)

        assert self.nmos is not None, ('printOrbitals(): '
                                   'No molecular orbitals have been collected')

        # Convert pol number to index
        o1 -= 1
        # Default to all.
        if o2 is None: o2 = self.nmos

        if space: print()

        # Print off orbitals in order of energy
        for n in range(o1, o2):

            # If this mo is not on a list, skip
            if MOs:
                try:
                    MOs.index(self.orbital_ids[n])
                except ValueError:
                    continue

            # The header for this MO
            head = 'Orbital {0}, Energy {1:7g} {2}, Occupation {3}'.format(
                    self.orbital_ids[n],
                    conv[units] * self.orbital_energies[n],
                    units,
                    self.orbital_occ[n])
            # The length of the lines above and below header
            print('    {0}'.format('-'*len(head)), head, '-'*len(head),
                                                                  sep='\n    ')
            # Table info header
            print('            Atomic Orb.  Symmetry    Percent')
            # Set the AO contubution sum to zero
            aosum = 0
            # Run over each AO to print its information
            aos = self.atomic_orbitals[n]
            for i in range(len(aos)):
                print('{0:12}{1}{0:4}{2}{0:2}{3:7.2%}'.format(' ',
                   aos.ao_id[i].center(8), aos.sym[i].center(10), aos.pcent[i]))
                # If the sum of all AOs is greater than the threshold, finish
                aosum += aos.pcent[i];
                if (aosum * 100 > thres) and thres != 1:
                    break

            # If we didn't go to 100%, tell what percent these AO make up of the MO
            if thres != 1 or self._ct :
                print('{0:11}{1}'.format(' ', '-'*33))
                print('{0:29}Total :{1:7.2%}'.format(' ', aosum))

            # One last seperator
            print('{0:4}{1}\n'.format(' ', '-'*len(head)))


    def printExcitations(self, order='energy', mod='all', trans=False,
                         thres=80, orb=False, units='eV'):
        '''Pretty-print the :py:attr:`excitation energies <excitation_energies>` 
        to the standard output.

        The excitation energies are printed in the *order* that you choose.
        The ordering schemes are:
        
         - 'energy' from lowest to highest :py:attr:`excitation_energies`
         - 'strength' from lowest to highest :py:attr:`oscillator_strengths`
         - 'number' from lowest to highest excitation number
         
        *mod* and can be:
         
         - the number of excitations to print
         - 'all' for all excitations
         - a symmetry group to print
        
        If the boolian *trans* is switched on, then after each
        excitation the :py:attr:`transitions` in that excitation will will be
        printed.  The number of transitions printed off is determined
        by *thres*, which is the maximum sum of the percent
        contributions to that excitation from the transitions that you
        want to be printed, starting with the transition with the
        largest percent contribution.  The boolian *orb* will return 
        the :py:func:`orbitals <printOrbitals>` that appear in the printed
        transitions.

        Raises :py:exc:`AssertionError` if no excitations are present.
        
        '''
        from .constants import HART2NM, HART2EV
        from numpy import argsort
        import re

        assert 'EXCITATIONS' in self.calctype, ('printExcitations(): '
                                         'No excitations have been collected.')
                                         
        if not re.search(r'au|eV', units):
            raise ChemDataError ('Invalid choice for units: ' + units)
        
        conv = {'au': 1, 'eV':HART2EV}

        # Determine what the title will say and how to set up the index list
        # based on the input parameters.
        #
        # For energy, the index is simply in acending order, since the
        # excitations are already ordered by energy.  The index is truncated
        # if only printing a subset of the excitations.  For excitations of
        # of a symmetry, the index of the excitations of that symmetry are
        # saved.
        #
        # The same methods are used for strength, except that the index is
        # sorted against the oscillator strengths.  When sorting for symmetry
        # care must be taken not to lose the indexes of the chosen symmetry
        # group, since sorting the index array against the oscillator
        # strengths sill return a second index array which is a sorted version
        # of the first index array, not a sorted index array for the
        # excitations.
        #
        # If looking for a specific excitation, simply make the index array
        # contain only that index.

        if order == 'energy':
            end = 'ordered by energy'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.excitation_energies[
                                self._excite_index])
            elif isinstance(mod, int):
                # Default to nexcite if mod is greater
                if mod > self.nexcite: mod = self.nexcite

                title = str(mod) + ' excitations, ' + end
                index = argsort(self.excitation_energies[
                                self._excite_index])
                index = index[:mod]
            else:
                title = 'All excitations of symmetry ' + mod + ', ' + end
                index = [i for i in self._excite_index 
                                       if self.excitation_symmetries[i] == mod]
                if not index:
                    raise ChemDataError (
                                       'Symmetry group ' + mod + ' not found!')
        elif order == 'strength':
            end = 'ordered by strength'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.oscillator_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
            elif isinstance(mod, int):
                # Default to nexcite if mod is greater
                if mod > self.nexcite: mod = self.nexcite

                title = str(mod) + ' excitations, ' + end
                index = argsort(self.oscillator_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
                index = index[:mod]
            else:
                title = 'All excitations of symmetry ' + mod + ', ' + end
                index = [i for i in self._excite_index
                                       if self.excitation_symmetries[i] == mod]
                if not index:
                    raise ChemDataError (
                                       'Symmetry group ' + mod + ' not found!')
                # Sort the index by strength, retaining the original numbers
                index = [index[i] for i in argsort(
                                             self.oscillator_strengths[index])]
                index = index[::-1] # reverse order
        elif order == 'number':
            title = 'Excitation number ' + str(mod)
            index = [int(mod)-1]
        else:
            raise ChemDataError ("Unknown keyword for 'order': " + order)

        # Title.  The title has ='s the length of the string on top and bottom
        title = title.join(('\n# ', ' #\n'))
        print('\n', title.center(3 * len(title) - 2, '='), sep='', end='\n\n')

        # Define heading
        if units == 'au':
            heading = ('\nEnergy (au)   Osc. Strength   Ex. Number   '
                   'Sym. Group   Wavelength (nm)\n')            
        else:
            heading = ('\nEnergy (ev)   Osc. Strength   Ex. Number   '
                   'Sym. Group   Wavelength (nm)\n')

        # Print once if not printing transitions
        if not trans:
            print(heading.center(71 * 2 + len(heading), '-'))

        # Print off the energies.  
        # Print the heading for each excitation if transitions are included.
        for i in self._excite_index[index]:

            if trans:
                print(heading.center(71 * 2 + len(heading), '-'))
            
            

            print('{0:9.5f}{5:4}{1:13.4E}{5:4}{2:>6}{5:7}{3}{5:6}{4:7.2f}'
                  .format((conv[units] * self.excitation_energies[i]),
                          self.oscillator_strengths[i],
                          str(i+1),
                          self.excitation_symmetries[i].center(10),
                          (HART2NM / self.excitation_energies[i]), ''))

            # Print off the transitions for this excitation if requested.
            # Spacing is determined similar to above.
            if trans:
                if orb:
                    orbs = []
                transsum = 0
                print('-' * 71)
                print('     Occ. MO       Unocc. MO     %Weight')
                print('   ', '-' * 36)
                for j in range(len(self.transitions[i])):
                    # Print heading for this excitation
                    occ = self.transitions[i].occ[j]
                    unocc = self.transitions[i].unocc[j]
                    comp = self.transitions[i].pcent[j]
                    # Add the MOs to the orb list unless it already exists
                    if orb and not orbs.count(occ):
                        orbs.append(occ)
                    if orb and not orbs.count(unocc):
                        orbs.append(unocc)
                    print('     {0}{1}{2}     {3:5.2%}'.format(occ.center(8),
                                            ' ->   ', unocc.center(10), comp))
                    # Sum up the percent composition to compare to threshold
                    transsum += comp
                    if (transsum * 100 > thres) and thres != 1:
                        break

                if thres != 1 or self._ct:
                    print('   ', '-' * 36)
                    print(' ' * 25, 'Total : {0:5.2%}'.format(transsum))
                    if not orb:
                        print('   ', '-' * 36)
                        print()
                    else:
                        self.printOrbitals(MOs=orbs, thres=thres, space=False)

        if not trans:
            print('-' * 71, end='\n\n')

    def printCircularDichroism(self, order='energy', mod=10, units='eV'):
        '''Pretty-print the circular dichroism results to the standard output.

        The excitation energies are printed in the *order* that you choose.
        The ordering schemes are:
        
         - 'energy' from lowest to highest :py:attr:`excitation_energies`
         - 'strength' from lowest to highest :py:attr:`oscillator_strengths`
         - 'number' from lowest to highest excitation number
         
        *mod* and can be:
         
         - the number of excitations to print
         - 'all' for all excitations
         - a symmetry group to print
        
        Raises :py:exc:`AssertionError` if no circular dichroism results are present.
        
        '''
        from .constants import HART2NM, HART2EV
        from numpy import argsort
        import re

        assert 'CD SPECTRUM' in self.calctype, ('printCircularDichroism(): '
                                         'No circular dichroism results have been collected.')
                                         
        if not re.search(r'au|eV', units):
            raise ChemDataError ('Invalid choice for units: ' + units)
        
        conv = {'au': 1, 'eV':HART2EV}

        if order == 'energy':
            end = 'ordered by energy'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.excitation_energies[
                                self._excite_index])
            elif isinstance(mod, int):
                # Default to nexcite if mod is greater
                if mod > self.nexcite: mod = self.nexcite

                title = str(mod) + ' excitations, ' + end
                index = argsort(self.excitation_energies[
                                self._excite_index])
                index = index[:mod]
            else:
                title = 'All excitations of symmetry ' + mod + ', ' + end
                index = [i for i in self._excite_index 
                                       if self.excitation_symmetries[i] == mod]
                if not index:
                    raise ChemDataError (
                                       'Symmetry group ' + mod + ' not found!')
        elif order == 'strength':
            end = 'ordered by strength'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.oscillator_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
            elif isinstance(mod, int):
                # Default to nexcite if mod is greater
                if mod > self.nexcite: mod = self.nexcite

                title = str(mod) + ' excitations, ' + end
                index = argsort(self.oscillator_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
                index = index[:mod]
            else:
                title = 'All excitations of symmetry ' + mod + ', ' + end
                index = [i for i in self._excite_index
                                       if self.excitation_symmetries[i] == mod]
                if not index:
                    raise ChemDataError (
                                       'Symmetry group ' + mod + ' not found!')
                # Sort the index by strength, retaining the original numbers
                index = [index[i] for i in argsort(
                                             self.oscillator_strengths[index])]
                index = index[::-1] # reverse order
        elif order == 'number':
            title = 'Excitation number ' + str(mod)
            index = [int(mod)-1]
        else:
            raise ChemDataError ("Unknown keyword for 'order': " + order)

        # Title.  The title has ='s the length of the string on top and bottom
        title = title.join(('\n# ', ' #\n'))
        print('\n', title.center(3 * len(title) - 2, '='), sep='', end='\n\n')
        print('Rotatory strengths are given in 1E-40 esu^2 cm^2', end='\n\n')

        # Define heading
        if units == 'au':
            heading = ('\nEnergy (au)   Rotatory Strength   Ex. Number   '
                   'Sym. Group   Wavelength (nm)\n')            
        else:
            heading = ('\nEnergy (ev)   Rotatory Strength   Ex. Number   '
                   'Sym. Group   Wavelength (nm)\n')

        # Print once if not printing transitions
        print(heading.center(75 * 2 + len(heading), '-'))

        # Print off the energies.  
        for i in self._excite_index[index]:

            print('{0:9.5f}{5:6}{1:13.4E}{5:6}{2:>6}{5:7}{3}{5:6}{4:7.2f}'
                  .format((conv[units] * self.excitation_energies[i]),
                          self.opt_rot_strengths[i],
                          str(i+1),
                          self.excitation_symmetries[i].center(10),
                          (HART2NM / self.excitation_energies[i]), ''))

        print('-' * 75, end='\n\n')

    def printTPA(self, order='energy', mod='all', units='eV'):
        '''Pretty-print the two-photon absorbance data to the standard 
        output.  

        The excitation energies and two-photon absorbance data are printed 
        in the *order* that you choose.  The ordering schemes are:
        
         - 'energy' from lowest to highest :py:attr:`excitation_energies`
         - 'linstrength' from lowest to highest :py:attr:`linear_tpa_strengths`
         - 'cirstrength' from lowest to highest :py:attr:`circular_tpa_strengths`
         - 'number' from lowest to highest excitation number
         
        *mod* and can be:
         
         - the number of excitations to print
         - 'all' for all excitations
        
        Raises :py:exc:`AssertionError` if no two-photon absorbance data is present.
        
        '''
        from .constants import HART2EV
        from numpy import argsort
        import re

        assert 'TPA' in self.calctype, ('printTPA(): '
                                        'No two-photon absorbance data '
                                        'has been collected.')

        if not re.search(r'au|eV', units):
            raise ChemDataError ('Invalid choice for units: ' + units)
        
        conv = {'au': 1, 'eV':HART2EV}
  
        if order == 'energy':
            end = 'ordered by energy'
            if mod == 'all':
                title = 'All excitations, ' + end 
                index = argsort(self.excitation_energies[
                                self._excite_index])
        elif order == 'linstrength':
            end = 'ordered by linear TPA strength'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.linear_tpa_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
        elif order == 'cirstrength':
            end = 'ordered by circular TPA strength'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.circular_tpa_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
        elif order == 'number':
            title = 'Excitation number ' + str(mod)
            index = [int(mod)-1]
        else:
            raise ChemDataError ("Unknown keyword for 'order': " + order)

        # Title.  The title has ='s the length of the string on top and bottom
        title = title.join(('\n# ', ' #\n'))
        print('\n', title.center(3 * len(title) - 2, '='), sep='', end='\n\n')

        print('TPA strengths are given in a.u.', end='\n\n')

        # Define heading
        if units == 'au':
            heading = ('\nEnergy (au)   Lin. TPA Strength   Cir. TPA Strength   '
                   'Ex. Number   Sym. Group\n')
        else:
            heading = ('\nEnergy (ev)   Lin. TPA Strength   Cir. TPA Strength   '
                   'Ex. Number   Sym. Group\n')
                 
        print(heading.center(77 * 2 + len(heading), '-'))

        # Print off the energies and two-photon absorbance strengths.  
        for i in self._excite_index[index]:

            print('{0:9.5f}{5:5}{1:13.4E}{5:7}{2:13.4E}{5:6}{3:>6}{5:8}{4}'
                  .format((conv[units] * self.excitation_energies[i]),
                          self.linear_tpa_strengths[i],
                          self.circular_tpa_strengths[i],
                          str(i+1),
                          self.excitation_symmetries[i].center(10),
                          ''))

        print('-' * 77, end='\n\n')

    def print3PA(self, order='energy', mod='all', units='eV'):
        '''Pretty-print the three-photon absorbance data to the standard 
        output.  

        The excitation energies and two-photon absorbance data are printed 
        in the *order* that you choose.  The ordering schemes are:
        
         - 'energy' from lowest to highest :py:attr:`excitation_energies`
         - 'linstrength' from lowest to highest :py:attr:`linear_3pa_strengths`
         - 'cirstrength' from lowest to highest :py:attr:`circular_3pa_strengths`
         - 'number' from lowest to highest excitation number
         
        *mod* and can be:
         
         - the number of excitations to print
         - 'all' for all excitations
        
        Raises :py:exc:`AssertionError` if no three-photon absorbance data is present.
        
        '''
        from .constants import HART2EV
        from numpy import argsort
        import re

        assert '3PA' in self.calctype, ('print3PA(): '
                                        'No three-photon absorbance data '
                                        'has been collected.')

        if not re.search(r'au|eV', units):
            raise ChemDataError ('Invalid choice for units: ' + units)
        
        conv = {'au': 1, 'eV':HART2EV}
  
        if order == 'energy':
            end = 'ordered by energy'
            if mod == 'all':
                title = 'All excitations, ' + end 
                index = argsort(self.excitation_energies[
                                self._excite_index])
        elif order == 'linstrength':
            end = 'ordered by linear 3PA strength'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.linear_3pa_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
        elif order == 'cirstrength':
            end = 'ordered by circular 3PA strength'
            if mod == 'all':
                title = 'All excitations, ' + end
                index = argsort(self.circular_3pa_strengths[
                                self._excite_index])
                index = index[::-1] # reverse order
        elif order == 'number':
            title = 'Excitation number ' + str(mod)
            index = [int(mod)-1]
        else:
            raise ChemDataError ("Unknown keyword for 'order': " + order)

        # Title.  The title has ='s the length of the string on top and bottom
        title = title.join(('\n# ', ' #\n'))
        print('\n', title.center(3 * len(title) - 2, '='), sep='', end='\n\n')

        print('3PA strengths are given in a.u.', end='\n\n')

        # Define heading
        if units == 'au':
            heading = ('\nEnergy (au)   Lin. 3PA Strength   Cir. 3PA Strength   '
                   'Ex. Number   Sym. Group\n')
        else:
            heading = ('\nEnergy (ev)   Lin. 3PA Strength   Cir. 3PA Strength   '
                   'Ex. Number   Sym. Group\n')
                 
        print(heading.center(77 * 2 + len(heading), '-'))

        # Print off the energies and two-photon absorbance strengths.  
        for i in self._excite_index[index]:

            print('{0:9.5f}{5:5}{1:13.4E}{5:7}{2:13.4E}{5:6}{3:>6}{5:8}{4}'
                  .format((conv[units] * self.excitation_energies[i]),
                          self.linear_3pa_strengths[i],
                          self.circular_3pa_strengths[i],
                          str(i+1),
                          self.excitation_symmetries[i].center(10),
                          ''))

        print('-' * 77, end='\n\n')
 
    def make_fragment_file(self, fragfile):
        '''Splits the system into user defined fragments and writes a fragment
        file for fragment analysis.
           
        '''
        frag_ele = []
        atoms = self.elements
        count = 0
        while len(atoms) > 0:
            s = ''
            for a in atoms:
                s += a + ' '
            print(s)
            inp = raw_input("Enter elements to put in fragment " + str(count+1) + ': ')
            inp = set(inp.split())
            if inp.issubset(atoms):
                frag_ele.append(inp)
                atoms.difference_update(inp)
                count += 1
            else:
                print("Not a valid set of elements")
        #print(frag_ele)
        
        from numpy import array
        
        frags = []
        for j in range(count):
            frags.append([])
        #print(frags)
        for i in range(self.natoms):
            for j in range(len(frag_ele)):
                if self.atoms[i] in frag_ele[j]:
                    frags[j].append('' + str(i+1) + ' ' + self.atoms[i])
        #print(frags)
        
        with open(fragfile + '.frag', 'w') as fl:
            for i in range(len(frags)):
                print('Fragment ' + str(i+1), file=fl)
                for j in frags[i]:
                    print(j, file=fl)


    def read_fragment_file(self, fragfile):
        '''Reads in the fragment file and returns a nested list of fragments
        This data is used for performing a fragment analysis in conjunction
        with the orbitals and excitations.

        It should be formatted like this:  
              
         | Fragment 1        
         | 1 C        
         | 2 C        
         | 3 C        
         | Fragment 2        
         | 4 C        
         | 5 C        
         | 6 C        

        Comments beginning with # are allowed.
        
        '''
        import re
        frags = []
        i = -1
        with open(fragfile) as fl:
            for line in fl:
                if line[0] == '#': # Skip a comment
                    continue
                # Verify the file was set up properly and increment counter
                if re.search('fragment', line, re.IGNORECASE):
                    if re.search(r'fragment \d', line, re.IGNORECASE):
                        frags.append([])
                        i += 1
                        continue
                    else:
                        raise ChemDataError('Fragment file set up incorrectly')

                # Store the atom or fragment
                if re.search(r'\d{1,3} \w+', line):
                    try:
                        frags[i].append(line.strip())
                    except IndexError:
                        raise ChemDataError('Fragment file set up incorrectly')

        return frags

    def fragment_analysis(self, frags, ex_type='ct', fragment=0,
                                        priority='lower', mthres=60, ethres=60):
        '''Method to determine which excitations belong to a certain
        excitation type.

        *frags* is a nested list of which atomic orbitals belong to which
        fragment.  It is the return value of :py:func:`read_fragment_file`.

        The choices for *ex_type* are:
        
         - 'localized'
         - 'ct'
         
        The argument *fragment* is required for type 'localized', and
        tells which fragment to localize on using the fragment number
        (1, 2, 3, etc...).  The number 0 means delocalized (i.e. is
        localized on zero fragments).

        *mthres* determines the threshold for an MO to be considered part of a
        fragment.
        
        *thres* determines the threshold for an excitation to be considered 
        part of a fragment.

        If no molecular orbitals or excitations have been collected, 
        :py:exc:`AssertionError` will be raised.
        
        '''
        from numpy import delete, s_, array
        import re

        assert self.nmos is not None, ('fragment_analysis(): '
                                   'No molecular orbitals have been collected')
        assert 'EXCITATIONS' in self.calctype, ('fragment_analysis(): '
                                              'No excitations were collected.')

        if not re.search(r'ct|localized', ex_type):
            raise ChemDataError ('Invalid choice for ex_type: ' + ex_type)

        if len(frags) < 2:
            raise ChemDataError ('Must define at least two fragments')
        if not re.search(r'higher|lower', priority):
            priority = 'lower'
            

        # Determine the fragments
        MO_frags = {}
        # Make the frags array into a dictionary telling the location of
        # each atom
        loc = dict([(j, i+1) for i in range(len(frags)) for j in frags[i]])

        # Make a list temporarily so it is mutable
        self.atomic_orbitals = list(self.atomic_orbitals)
        # Loop over each MO to determine which fragment it is localized on
        for n in range(self.nmos):
            fragsum = {}
            for i in range(len(frags)):
                fragsum[i+1] = 0 # Initiallize sums
            aos = self.atomic_orbitals[n]
            for i in range(len(aos)):
                # Add this contribution to the correct frag
                if aos.ao_id[i] in loc:
                    fragsum[loc[aos.ao_id[i]]] += aos.pcent[i]
                else:
                    raise ChemDataError ('Missing ' + str(aos.ao_id[i]) +
                                                     ' from the fragment file')
            # If the sum is over the threshold, then assign fragment
            for i in range(len(frags)):
                if fragsum[i+1]*100 > mthres:
                    MO_frags[self.orbital_ids[n]] = i+1
                    if priority == 'lower':
                        break

        # Make tuple once more to be immutable
        self.atomic_orbitals = tuple(self.atomic_orbitals)
        self._ct = True

        # Determine which excitations belong to chosen type
        self.transitions = list(self.transitions)
        remove = []
        # Loop over all excitations
        for n in self._excite_index:
            transsum = 0
            # Loop over each transition in the excitatiopn
            for i in range(len(self.transitions[n])):
                occ = self.transitions[n].occ[i]
                unocc = self.transitions[n].unocc[i]
                # If one of the orbitals did not make the cut at all, skip
                if occ not in MO_frags or unocc not in MO_frags:
                    continue
                if ex_type == 'ct':
                    if MO_frags[occ] == MO_frags[unocc]:
                        # If the fragments are different the same, it's not CT.
                        continue
                    if not (MO_frags[occ] or MO_frags[unocc]):
                        # If they are both 0, it is not CT
                        continue
                elif ex_type == 'localized':
                    # If the fragments are different, it's not localized
                    if MO_frags[occ] != MO_frags[unocc]:
                        continue
                    if fragment:
                        # If one of fragments isn't what we want, chuck it 
                        # The above test ensures both fragments are the same
                        if MO_frags[occ] != fragment:
                            continue
                # We've passed the tests, now we'll sum up the contributions
                transsum += self.transitions[n].pcent[i]
                # If we make it above the cutoff then this
                # excitation passes the test. 
                if transsum*100 > ethres:
                    # Truncate the transitions list to applicable transitions
                    self.transitions[n] = delete(self.transitions[n],
                                           s_[i+1:len(self.transitions[n])])
                    break
            else:
                # This excitation did not pass the tests
                remove.append(n)

        # Remove excitations that are not of the type we are interested in
        remove.reverse()
        from numpy import delete
        for i in remove:
            self._excite_index = delete(self._excite_index, i)
        self.transitions = tuple(self.transitions)

        if not len(self._excite_index):
            raise ChemDataError ('No excitations matching your '
                                                        'criteria were found.')

    def sum_over_states(self, maxval=None, gamma=None, omega=None):
        '''Uses the sum over states formula to determine the polarizability'''
        from numpy import array, asarray, zeros

        assert 'EXCITATIONS' in self.calctype, ('sum_over_states():'
                                              'No excitations were collected.')
        if ((gamma is None and omega is not None) or 
            (gamma is not None and omega is None)):
            raise AssertionError ('sum_over_states: either both or neither of '
                                  'gamma and omega must be given')

        # Determine the number of frequencies from omega.  It is 1 if static
        # Make sure all values are a numpy array
        if omega is None:
            self.calctype.update(['POLARIZABILITY', 'STATIC'])
            self.npol = 1
            self.e_frequencies = array([0.0])
            gamma = array([0.0])
        else:
            self.calctype.update(['POLARIZABILITY', 'FD'])
            gamma = array([gamma])
            # Omega is either a range of frequencies or a single frequency
            try:
                self.npol = len(omega)
                self.e_frequencies = asarray(omega)
            except TypeError:
                self.npol = 1
                self.e_frequencies = array([omega])

        # If no max is given, use the number of excitations
        if maxval is None:
            maxval = self.nexcite
        elif maxval > self.nexcite:
            from sys import stderr
            msg = 'Given maxval {0} is greater than nexcite {1}'
            print(msg.format(maxval, self.nexcite), file=stderr)
            maxval = self.nexcite

        def polcomp(maxsum, a, b, TDM, excite_energies, omega, Gamma):
            r'''
            Calculate a polarizability tensor component using the
            sum-over-states formula

            \alpha_{\alpha\beta} = \sum^N_{n \neq 0}
              \frac{\mu_{n\alpha}\mu_{n\beta}}{\omega_{0,n} - \omega - i\Gamma}
             +\frac{\mu_{n\alpha}\mu_{n\beta}}{\omega_{0,n} + \omega + i\Gamma}
            '''
            mu1, mu2 = TDM[:maxsum,a], TDM[:maxsum,b]
            omega_0 = excite_energies[:maxsum]
            return sum(mu1 * mu2 / ( omega_0 - omega - 1j * Gamma )
                     + mu1 * mu2 / ( omega_0 + omega + 1j * Gamma )
                      )

        # Set the polarizability values
        ee = self.excitation_energies
        self.qm_pol = zeros((self.npol,3,3), dtype=complex)
        for i, om in enumerate(self.e_frequencies):
            for a in range(3):
                for b in range(3):
                    self.qm_pol[i,a,b] = polcomp(maxval, a, b,
                                                         self.TDM, ee,
                                                         om, gamma)

        # If static, only real part is returned
        if 'STATIC' in self.calctype:
            self.qm_pol = self.qm_pol.real

        # Return the number summed in case we want it
        return maxval

    def sum_over_states_gamma(self, maxval=None, gamma=None, omega=None):
        '''Uses the sum over states formula to determine the negative part of
           the second hyperpolarizability'''
        from numpy import array, asarray, zeros, einsum

        assert 'EXCITATIONS' in self.calctype, ('sum_over_states():'
                                              'No excitations were collected.')
        if ((gamma is None and omega is not None) or 
            (gamma is not None and omega is None)):
            raise AssertionError ('sum_over_states: either both or neither of '
                                  'gamma and omega must be given')

        # Determine the number of frequencies from omega.  It is 1 if static
        # Make sure all values are a numpy array
        if omega is None:
            self.calctype.update(['SECOND HYPERPOLARIZABILITY', 'STATIC'])
            self.npol = 1
            self.e_frequencies = array([0.0])
            gamma = 0.0
        else:
            self.calctype.update(['SECOND HYPERPOLARIZABILITY', 'FD'])
            gamma = gamma
            # Omega is either a range of frequencies or a single frequency
            try:
                self.npol = len(omega)
                self.e_frequencies = asarray(omega)
            except TypeError:
                self.npol = 1
                self.e_frequencies = array([omega])

        # If no max is given, use the number of excitations
        if maxval is None:
            maxval = self.nexcite
        elif maxval > self.nexcite:
            from sys import stderr
            msg = 'Given maxval {0} is greater than nexcite {1}'
            print(msg.format(maxval, self.nexcite), file=stderr)
            maxval = self.nexcite

        def shpolcomp(maxsum, a, b, c, d, TDM, excite_energies, omega, Gamma):
            r'''
            Calculate the pure one-photon part of a second hyperpolarizability
            tensor component using the sum-over-states formula

            \gamma_{abcd} = \sum^k_{i \neq 0}\sum^n_{j \neq 0}......
            '''
            mu1, mu2, mu3, mu4 = TDM[:maxsum,a], TDM[:maxsum,b], TDM[:maxsum,c], TDM[:maxsum,d]
            omega_0 = excite_energies[:maxsum]

            #a1 = 1.0 / (omega_0 - omega_b - omega_c - omega_d - 1j * Gamma)
            #a2 = 1.0 / (omega_0 - omega_b - omega_c           - 1j * Gamma)
            #a3 = 1.0 / (omega_0 - omega_b - omega_d           - 1j * Gamma)
            #a4 = 1.0 / (omega_0 - omega_c - omega_d           - 1j * Gamma)
            #a5 = 1.0 / (omega_0 - omega_b                     - 1j * Gamma)
            #a6 = 1.0 / (omega_0 - omega_c                     - 1j * Gamma)
            #a7 = 1.0 / (omega_0 - omega_d                     - 1j * Gamma)

            #b1 = 1.0 / (omega_0 + omega_b + omega_c + omega_d + 1j * Gamma)
            #b2 = 1.0 / (omega_0 + omega_b + omega_c           + 1j * Gamma)
            #b3 = 1.0 / (omega_0 + omega_b + omega_d           + 1j * Gamma)
            #b4 = 1.0 / (omega_0 + omega_c + omega_d           + 1j * Gamma)
            #b5 = 1.0 / (omega_0 + omega_b                     + 1j * Gamma)
            #b6 = 1.0 / (omega_0 + omega_c                     + 1j * Gamma)
            #b7 = 1.0 / (omega_0 + omega_d                     + 1j * Gamma)
            a1 = 1.0 / (omega_0 - 1*omega - 1j * Gamma)
            a2 = 1.0 / (omega_0 - 2*omega - 1j * Gamma)
            a3 = 1.0 / (omega_0 - 0*omega - 1j * Gamma)
            a4 = 1.0 / (omega_0 - 0*omega - 1j * Gamma)
            a5 = 1.0 / (omega_0 - 1*omega - 1j * Gamma)
            a6 = 1.0 / (omega_0 - 1*omega - 1j * Gamma)
            a7 = 1.0 / (omega_0 + 1*omega - 1j * Gamma)

            b1 = 1.0 / (omega_0 + 1*omega + 1j * Gamma)
            b2 = 1.0 / (omega_0 + 2*omega + 1j * Gamma)
            b3 = 1.0 / (omega_0 + 0*omega + 1j * Gamma)
            b4 = 1.0 / (omega_0 + 0*omega + 1j * Gamma)
            b5 = 1.0 / (omega_0 + 1*omega + 1j * Gamma)
            b6 = 1.0 / (omega_0 + 1*omega + 1j * Gamma)
            b7 = 1.0 / (omega_0 - 1*omega + 1j * Gamma)

            temp = 0.0
            temp -= einsum('i,i,j,j,i,i,j', mu1, mu2, mu3, mu4, a1, a5, a7)
            temp -= einsum('i,i,j,j,i,i,j', mu1, mu2, mu4, mu3, a1, a5, a6)
            temp -= einsum('i,i,j,j,i,i,j', mu1, mu3, mu2, mu4, a1, a6, a7)
            temp -= einsum('i,i,j,j,i,i,j', mu1, mu3, mu4, mu2, a1, a6, a5)
            temp -= einsum('i,i,j,j,i,i,j', mu1, mu4, mu2, mu3, a1, a7, a6)
            temp -= einsum('i,i,j,j,i,i,j', mu1, mu4, mu3, mu2, a1, a7, a5)

            temp -= einsum('i,i,j,j,i,j,j', mu1, mu2, mu3, mu4, a5, b6, a7)
            temp -= einsum('i,i,j,j,i,j,j', mu1, mu2, mu4, mu3, a5, b7, a6)
            temp -= einsum('i,i,j,j,i,j,j', mu1, mu3, mu2, mu4, a6, b5, a7)
            temp -= einsum('i,i,j,j,i,j,j', mu1, mu3, mu4, mu2, a6, b7, a5)
            temp -= einsum('i,i,j,j,i,j,j', mu1, mu4, mu2, mu3, a7, b5, a6)
            temp -= einsum('i,i,j,j,i,j,j', mu1, mu4, mu3, mu2, a7, b6, a5)

            temp -= einsum('i,i,j,j,i,i,j', mu2, mu1, mu3, mu4, b1, b5, b6)
            temp -= einsum('i,i,j,j,i,i,j', mu2, mu1, mu4, mu3, b1, b5, b7)
            temp -= einsum('i,i,j,j,i,i,j', mu3, mu1, mu2, mu4, b1, b6, b5)
            temp -= einsum('i,i,j,j,i,i,j', mu3, mu1, mu4, mu2, b1, b6, b7)
            temp -= einsum('i,i,j,j,i,i,j', mu4, mu1, mu2, mu3, b1, b7, b5)
            temp -= einsum('i,i,j,j,i,i,j', mu4, mu1, mu3, mu2, b1, b7, b6)

            temp -= einsum('i,i,j,j,i,j,j', mu2, mu1, mu3, mu4, b5, a7, b6)
            temp -= einsum('i,i,j,j,i,j,j', mu2, mu1, mu4, mu3, b5, a6, b7)
            temp -= einsum('i,i,j,j,i,j,j', mu3, mu1, mu2, mu4, b6, a7, b5)
            temp -= einsum('i,i,j,j,i,j,j', mu3, mu1, mu4, mu2, b6, a5, b7)
            temp -= einsum('i,i,j,j,i,j,j', mu4, mu1, mu2, mu3, b7, a6, b5)
            temp -= einsum('i,i,j,j,i,j,j', mu4, mu1, mu3, mu2, b7, a5, b6)

            return temp 

        # Set the polarizability values
        ee = self.excitation_energies
        self.secondhyperpolarizability = zeros((self.npol,3,3,3,3),dtype='complex128')
        for i, om in enumerate(self.e_frequencies):
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        for d in range(3):
                            self.secondhyperpolarizability[i,a,b,c,d] = shpolcomp(maxval, a, b, c, d,
                                                                                  self.TDM, ee,
                                                                                  om, gamma)
        # If static, only real part is returned
        if 'STATIC' in self.calctype:
            self.secondhyperpolarizability = self.secondhyperpolarizability.real

        # Return the number summed in case we want it
        return maxval
