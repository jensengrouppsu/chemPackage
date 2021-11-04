from __future__ import print_function, division
from .errorclass import ChemDataError
from numpy import array, append, reshape, bincount

class Coordinates(object):
    '''Extends ChemData class with methods to manipulate coordinates
    and normal modes.
    '''

    def printCoords(self, mode=None, dim=False, qm=True, a1=1, a2=None,
                    file=None, dalton=False, abasis=None, latex=False):
        '''Prints the geometry coordinates to screen.

        *mode* determines how the coordinates will be printed.  If
        omitted, the coordinates will be printed without numbers.

        Valid options are:
         - 'num':       Prints the number of each atom with the element.
         - 'xyz':       Prints the total number of atoms, then a space, then the coordinates.
         - 'xyz_title': Same as 'xyz', but prints the title instead of a space, if a title is available.
         - 'dimblock':  Same as 'xyz', but omits the space altogether.

        The *dim* and *qm* booleans control what set of coordinates
        are written. If both are true, the QM coordinates will be written
        before the DIM coordinates.

        *a1* is the lowest atom to print
        *a2* is the highest atom, defaulting to the number of atoms.
        This option allows you to print a subset of atoms.
        Note that if you include both QM and DIM atoms, the indexes
        reflect the two systems together.

        *file* allows you to select where to print to.  If omitted, it
        will print to standard error.  You may give the name of a file
        or an already open file object to write to.  Or, if file is
        'xyz', then it will create a .xyz file based on the filename
        attribute and print to that.

        *dalton* allows the coordinates to be printed in Dalton input
        format, with the atoms separated by element and given the correct
        headings.

        *abasis* allows the program to read in the atom basis for Dalton
        calculations where each element is represented by a different basis
        set.

        *latex* allows the output to be printed in a LaTeX table type of
        format.

        Raises and AssertionError if *dim* is True and not a
        DIMQM calculation, if an invalid mode is entered, or if file
        is 'xyz' and the filename attribute does not exist

        '''
        import os, sys
        print(mode, dim, qm, a1, a2, file, dalton, abasis, latex)

        if dim: assert 'DIM' in self.calctype, ('printCoords(): '
                                                    'Not a DIMQM calculation.')
        if file == 'xyz': assert 'filename' in self, ('printCoords(): '
                                              'filename attribute not defined')
        assert (mode is None or
                mode in ('xyz', 'xyz_title', 'num', 'dimblock')), (
                                          'printCoords(): Invalid mode: '+mode)

        # Create an xyz file if file is True
        if file == 'xyz':
            # Split current extention off filename, then replace with .xyz
            file = '.'.join([os.path.splitext(self.filename)[0], 'xyz'])
            # Open, and remember we need to close after
            file = open(file, 'w')
            closebool = True
        # If file is None, use standard out
        elif file is None:
            file = sys.stdout
            closebool = False
        # Otherwise, try to open the file
        else:
            try:
                file = open(file, 'w')
            # If it fails, then it was already an open file
            except TypeError:
                closebool = False
            # If it suceeds, then remmeber that we must close the file
            else:
                closebool = True

        # Prep the coordinates and atoms
        if dim and qm:
            coords = self.allcoords
            atoms = self.allatoms
        elif dim:
            coords = self.dim_coordinates
            atoms = self.dim_atoms
        else:
            coords = self.coordinates
            atoms = self.atoms

        # For Dalton formatted printing, we need to sort the atoms and
        # coordinates.
        if dalton:
            from .constants import elem
            sorted_atoms = []
            sorted_coords = array([], dtype=float)
            nuclear_charge = {}
            for item in self.elements:
                nuclear_charge[item] = str(float(elem.index(item)))
                for i in range(self.natoms):
                    if atoms[i] == item:
                        sorted_atoms = append(sorted_atoms, atoms[i])
                        sorted_coords = append(sorted_coords,
                                              [coords[i][0],coords[i][1],coords[i][2]])
            atoms = sorted_atoms
            coords = sorted_coords.reshape(-1,3)

        # Convert atom number to index
        a1 -= 1
        # Default to all.
        if a2 is None: a2 = len(atoms)

        # Determine how to format the line, and make numbers correct
        if mode == 'num':
            # Find string length of largest number to be printed and make that
            # number a string
            maxlen = str(len(str(( a2 - a1 ) + 1)))
            fmt = '{0:<'+maxlen+'} {1:<2}{2[0]:14.8f}{2[1]:14.8f}{2[2]:14.8f}'
            num = 1
        else:
            if latex == True:
                fmt = '{0:>3}{2:>4}{1[0]:12.6f}{2:>2}{1[1]:12.6f}{2:>2}{1[2]:12.6f}{3:<2}'
            else:
                fmt = '{0:<2}{1[0]:14.8f}{1[1]:14.8f}{1[2]:14.8f}'

        # Make number of atoms
        natoms = len(atoms[a1:a2])

        # Print the number of atoms if the mode calls for it
        if mode in ('xyz', 'xyz_title', 'dimblock'): print(natoms, file=file)
        # Print a space if 'xyz', or if 'xyz_title' and there is no title
        if mode == 'xyz' or (mode == 'xyz_title' and 'title' not in self):
            print(file=file)
        # Print the title if 'xyz_title' and there is a title
        elif mode == 'xyz_title':
            print(self.title, file=file)
        # Print the coordinates
        if dalton == False:
            for i in range(a1, a2):
                if mode == 'num':
                    print(fmt.format(num, atoms[i], coords[i]), file=file)
                    num += 1
                else:
                    if latex == True:
                        print(fmt.format(atoms[i], coords[i], '&', '\\\\'), file=file)
                        print(' \hline', file=file)
                    else:
                        print(fmt.format(atoms[i], coords[i]), file=file)
        # Print the coordinates in Dalton style
        elif dalton == True:
            for item in self.elements:
                lcharge = str(len(nuclear_charge[item]))
                lnum = str(list(atoms).count(item))
                lnum = str(len(lnum))
                if abasis is None:
                    fm = '{0}={1:<'+lcharge+'} {2}={3:<'+lnum+'}'
                    print(fm.format("Charge", nuclear_charge[item],
                                    "Atoms", list(atoms).count(item)), file=file)
                else:
                    lbf = str(len(abasis[item]))
                    fm = '{0}={1:<'+lcharge+'} {2}={3:<'+lnum+'} {4}={5:<'+lbf+'}'
                    print(fm.format("Charge", nuclear_charge[item],
                                    "Atoms", list(atoms).count(item),
                                    "Basis", abasis[item]), file=file)
                for i in range(a1, a2):
                    if atoms[i] == item:
                        print(fmt.format(atoms[i], coords[i]), file=file)

        # Close the file if appropriate
        if closebool: file.close()


    def writeCoords(self, dim=False, qm=True, a1=1, a2=None):
        '''This is a shortcut for :py:func:`printCoords(mode='xyz', file='xyz') <printCoords>`.

        This will write the coordinates to a .xyz file based on the
        name in the filename attribute.  Returns :py:exc:`AttributeError` if
        no filename exists.

        '''
        self.printCoords(mode='xyz', file='xyz', dim=dim, qm=qm, a1=a1, a2=a2)


    def writePDB(self, dim=False, qm=True, a1=1, a2=None):
        '''Writes the corrdinates to a .pdb file with the same name
        as the current file.

        The *dim* and *qm* booleans control what set of coordinates
        are written. If both are true, the QM coordinates will be written
        before the DIM coordinates.

        *a1* is the lowest atom to print
        *a2* is the highest atom, defaulting to the number of atoms.
        This option allows you to print a subset of atoms.
        Note that if you include both QM and DIM atoms, the indexes
        reflect the two systems together.
        '''
        import os
        from numpy import where, vstack

        # Define title format
        ft = 'CMPND  {0}'
        # Define coordinate format
        fc = ('HETATM{0:>5d}{1:>3}   LIG     1    '
              '{2[0]:>8.3f}{2[1]:>8.3f}{2[2]:>8.3f}  1.00  0.00          '
              '{1:>2}  ')
    
        # Var No.3 is the Hirshfeld induced charges (imag, z), 
        # filled into the temp. column. -- Pengchong Liu, 08/10/2016
        #fc = ('HETATM{0:>5d}{1:>3}   LIG     1    '
        #      '{2[0]:>8.3f}{2[1]:>8.3f}{2[2]:>8.3f}  1.00{3:>6.1f}          '
        #      '{1:>2}  ')
        # Grab Hirshfeld_induced_charges.imag at dir='Z'
        #hcharges = self.hirshfeld_induced_charges[:,2].imag
    
        # Define bond format
        fb = 'CONECT{0:>5d}'

        # Grab the correct coordinates
        if dim and qm:
            atoms = self.allatoms
            coords = self.allcoords
            bonds = self.bonds
            # Grab the DIM bonds, and add the number of QM atoms to each
            # element so that the indexes remain correct
            tmpbond = self.dim_bonds
            tmpbond += self.natoms
            # Append DIM bonds after QM bonds
            bonds = vstack((bonds, tmpbond))
            # Remember the total number of atoms in system
            totatoms = self.nallatoms
        elif dim:
            atoms = self.dim_atoms
            coords = self.dim_coordinates
            bonds = self.dim_bonds
            # Remember the total number of atoms in system
            totatoms = self.ndim
        else:
            atoms = self.atoms
            coords = self.coordinates
            bonds = self.bonds
            # Remember the total number of atoms in system
            totatoms = self.natoms
        
        # Convert atom number to index
        a1 -= 1
        # Default to all.
        if a2 is None: a2 = len(atoms)

        # Open the file as pdb with same name
        filename = '.'.join([os.path.splitext(self.filename)[0], 'pdb'])
        with open(filename, 'w') as fl:

            # pbd Starts with some info.  Give the title if there is one
            if self.title:
                title = ft.format('MOLECULE: '+self.title)
            else:
                title = ft.format('UNNAMED')
            print(title, file=fl)

            # Place the coordinates in the file
            i = 1
            #for atom, coord, hc in zip(atoms[a1:a2], coords[a1:a2], hcharges[a1:a2]):
            for atom, coord in zip(atoms[a1:a2], coords[a1:a2]):
                print(fc.format(i, atom, coord), file=fl)
                #print(fc.format(i, atom, coord, hc), file=fl)
                i += 1

            # Place the bonds in the file.  The bonding is redundant, so each
            # atom must be specified and bonds will be listed mutlitple times.
            for i in range(totatoms):

                # Skip atoms we don't want to show a bond to
                if i < a1 or i > a2: continue
                # Print the bonding keyword
                print(fb.format(i+1), end='', file=fl)
                # Find everywhere that this atoms is in the first column
                indx = where(bonds[:,0] == i)[0]
                # Print off all atoms this one is bonded to
                for j in indx:
                    # Skip atoms we don't want to show a bond to
                    if bonds[j,1] < a1 or bonds[j,1] > a2: continue
                    print('{0:>5d}'.format(bonds[j,1]+1), end='', file=fl)
                # Repeat for the second column
                indx = where(bonds[:,1] == i)[0]
                for j in indx:
                    # Skip atoms we don't want to show a bond to
                    if bonds[j,0] < a1 or bonds[j,0] > a2: continue
                    print('{0:>5d}'.format(bonds[j,0]+1), end='', file=fl)
                # New line
                print(file=fl)

            # End the file
            print('END   ', file=fl)

    def printGSGradients(self, a1=1, a2=None):
        '''Prints the ground state gradients to screen.

        *a1* is the lowest atom to print
        *a2* is the highest atom, defaulting to the number of atoms.
        This option allows you to print a subset of atoms.

        *file* allows you to select where to print to.  If omitted, it
        will print to standard error.  You may give the name of a file
        or an already open file object to write to.  Or, if file is
        'xyz', then it will create a .xyz file based on the filename
        attribute and print to that.

        '''
        import os, sys

        # Prep the gradients and atoms
        grad_types = self.gs_gradient.keys()
        grad = self.gs_gradient
        atoms = self.atoms

        # Convert atom number to index
        a1 -= 1
        # Default to all.
        if a2 is None: a2 = len(atoms)

        # Format the line for outputting gradients
        fmt = '{0:<2}{1[0]:12.6f}{1[1]:12.6f}{1[2]:12.6f}'

        # Ordering
        order = ('nuclear repulsion gradient',
                 'weighted density gradient',
                 'kinetic energy gradient',
                 '2-electron gradient',
                 'DFT CD+XC gradient',
                 'total DFT gradient',)

        # Table of headings
        table = {'nuclear repulsion gradient' : 'Nuclear Repulsion Gradient',
                 'weighted density gradient'  : 'Weighted Density Gradient',
                 'kinetic energy gradient'    : 'Kinetic Energy Gradient',
                 '2-electron gradient'        : 'Two-Electron Gradient',
                 'DFT CD+XC gradient'         : 'DFT Charge Density Fit + XC-Functional Gradient' ,
                 'total DFT gradient'         : 'Total DFT Ground State Gradient',}

        # Print the gradient
        print()
        print('~' * 42)
        print('All gradients are in Hartrees/Bohr (Eh/a0)')
        print('~' * 42)
        print()
        for item in order:
            if item in grad_types:
                print(table[item])
                print('=' * 47)
                for i in range(a1, a2):
                    print(fmt.format(atoms[i], grad[item][i]))
                print()

    def printESGradients(self, a1=1, a2=None, latex=False):
        '''Prints the excited state gradients to screen.

        *a1* is the lowest atom to print
        *a2* is the highest atom, defaulting to the number of atoms.
        This option allows you to print a subset of atoms.

        *file* allows you to select where to print to.  If omitted, it
        will print to standard error.  You may give the name of a file
        or an already open file object to write to.  Or, if file is
        'xyz', then it will create a .xyz file based on the filename
        attribute and print to that.

        *latex* allows you to output the data conveniently for putting
        it in a LaTeX data table.

        '''
        import os, sys

        # Prep the gradients and atoms
        grad_types = self.es_gradient.keys()
        grad = self.es_gradient
        atoms = self.atoms

        # Convert atom number to index
        a1 -= 1
        # Default to all.
        if a2 is None: a2 = len(atoms)

        # Format the line for outputting gradients
        if latex == False:
            fmt = '{0:<2}{1[0]:12.6f}{1[1]:12.6f}{1[2]:12.6f}'
        else:
            fmt = '{0:>3}{2:>4}{1[0]:12.6f}{2:>2}{1[1]:12.6f}{2:>2}{1[2]:12.6f}{3:<2}'

        # Ordering
        order = ('nuclear repulsion gradient',
                 'weighted density gradient',
                 'kinetic energy gradient',
                 '2-electron gradient',
                 'Vxc gradient',
                 'CD gradient',
                 'fxc gradient',
                 'TDDFT CD+XC gradient',
                 'TDDFT excitation energy gradient',
                 'total TDDFT gradient',)

        # Table of headings
        table = {'nuclear repulsion gradient'       : 'Nuclear Repulsion Gradient',
                 'weighted density gradient'        : 'Weighted Density Gradient',
                 'kinetic energy gradient'          : 'Kinetic Energy Gradient',
                 '2-electron gradient'              : 'Two-Electron Gradient',
                 'Vxc gradient'                     : 'XC-Potential Gradient',
                 'CD gradient'                      : 'Charge Density Fit Gradient',
                 'fxc gradient'                     : 'XC-Kernel Gradient',
                 'TDDFT CD+XC gradient'             : 'TDDFT Charge Density Fit + XC-Functional Gradient' ,
                 'TDDFT excitation energy gradient' : 'TDDFT Excitation Energy Gradient' ,
                 'total TDDFT gradient'             : 'Total TDDFT Excited State Gradient',}

        # Print the gradient
        print()
        print('~' * 42)
        print('All gradients are in Hartrees/Bohr (Eh/a0)')
        print('~' * 42)
        print()
        for item in order:
            if item in grad_types:
                print(table[item])
                print('=' * 47)
                for i in range(a1, a2):
                    if latex == False:
                        print(fmt.format(atoms[i], grad[item][i]))
                    else:
                        print(fmt.format(atoms[i], grad[item][i],'&','\\\\'))
                        print(' \hline')
                print()

    def __add__(self, other):
        '''Concatenates two molecules together and returns the new molecule.

        The attributes that are concatenated are:

          - :py:attr:`coordinates`
          - :py:attr:`atoms`
          - :py:attr:`natoms`
          - :py:attr:`dim_coordinates`
          - :py:attr:`dim_atoms`
          - :py:attr:`ndim`

        All other properties would be unphysical to concatenate
        and are emptied.

        '''
        from numpy import concatenate
        from .chemdata import ChemData
        new = ChemData()
        new.coordinates = concatenate((self.coordinates, other.coordinates))
        new.atoms = concatenate((self.atoms, other.atoms))
        new.natoms = len(new.atoms)
        if 'DIM' in self.calctype & other.calctype:
            new.dim_coordinates = concatenate((self.dim_coordinates,
                                               other.dim_coordinates))
            new.dim_atoms = concatenate((self.dim_atoms, other.dim_atoms))
            new.ndim = len(new.dim_atoms)
        elif 'DIM' in other.calctype:
            new.dim_coordinates = other.dim_coordinates.copy()
            new.dim_atoms = other.dim_atoms.copy()
            new.ndim = len(new.dim_atoms)
        elif 'DIM' in self.calctype:
            new.dim_coordinates = self.dim_coordinates.copy()
            new.dim_atoms = self.dim_atoms.copy()
            new.ndim = len(new.dim_atoms)

        return new


    def join(self, other):
        '''Concatenates one molecule into the current one.

        The attributes that are concatenated are:

          - :py:attr:`~.coordinates`
          - :py:attr:`~.atoms`
          - :py:attr:`~.natoms`
          - :py:attr:`~.elements`
          - :py:attr:`~.nelements`
          - :py:attr:`~.dim_coordinates`
          - :py:attr:`~.dim_atoms`
          - :py:attr:`~.ndim`
          - :py:attr:`~.dim_elements`
          - :py:attr:`~.ndim_elements`

        All other properties would be unphysical to concatenate
        and are emptied.

        '''
        from numpy import concatenate
        self.coordinates = concatenate((self.coordinates, other.coordinates))
        self.atoms = concatenate((self.atoms, other.atoms))
        self.natoms = len(self.atoms)
        self.elements.update(other.elements)
        self.nelements = len(self.elements)
        if 'DIM' in self.calctype & other.calctype:
            self.dim_coordinates = concatenate((self.dim_coordinates,
                                                other.dim_coordinates))
            self.dim_atoms = concatenate((self.dim_atoms, other.dim_atoms))
            self.ndim = len(self.dim_atoms)
            self.dim_elements.update(other.dim_elements)
            self.ndim_elements = len(self.dim_elements)
        elif 'DIM' in other.calctype:
            self.dim_coordinates = other.dim_coordinates.copy()
            self.dim_atoms = other.dim_atoms.copy()
            self.ndim = len(self.dim_atoms)
            self.dim_elements = other.dim_elements.copy()
            self.ndim_elements = len(self.dim_elements)

        self.empty(ignore=['coordinates', 'atoms', 'natoms',
                           'dim_coordinates', 'dim_atoms', 'ndim',
                           'elements', 'nelements', 'dim_elements',
                           'ndim_elements' ])

    # Remember that the property decorator makes a
    # method appear as an attribute
    @property
    def allcoords(self):
        '''This method is guarunteed to return the coordinates of all
        atoms in the system, QM or otherwise, with the QM atoms listed
        first.

        '''
        if 'DIM' in self.calctype:
            if self.coordinates is not None:
                from numpy import concatenate
                return concatenate((self.coordinates, self.dim_coordinates))
            else:
                return self.dim_coordinates
        else:
            return self.coordinates


    @property
    def allatoms(self):
        '''This method is guarunteed to return all atoms in the system,
        QM or otherwise, with the QM atoms listed first.

        '''
        if 'DIM' in self.calctype:
            if self.atoms is not None:
                from numpy import concatenate
                return concatenate((self.atoms, self.dim_atoms))
            else:
                return self.dim_atoms
        else:
            return self.atoms

    @property
    def allelements(self):
        '''This method is guarunteed to return elements in the system,
        QM or otherwise.

        '''
        if 'DIM' in self.calctype:
            if self.elements:
                return self.elements.union(self.dim_elements)
            else:
                return self.dim_elements
        else:
            return self.elements


    @property
    def nallatoms(self):
        '''This method is guarunteed to return the total atoms in the system,
        QM or otherwise.

        '''

        if 'DIM' in self.calctype:
            return self.natoms + self.ndim
        else:
            return self.natoms


    @property
    def nallelements(self):
        '''This method is guarunteed to return the total number of elements
        in the system, QM or otherwise.

        '''
        if 'DIM' in self.calctype:
            return len(self.elements + self.dim_elements)
        else:
            return self.nelements


    @property
    def masses(self):
        '''Returns the mass of each atom in QM system.'''
        from numpy import vectorize
        from .constants import atomic_mass

        # Make the atomic mass function broadcastable on a numpy array
        atm_mass = vectorize(atomic_mass)
        # Find mass for each atom
        return atm_mass(self.atoms)


    @property
    def atomic_numbers(self):
        '''Returns the atomic number of each atom in the QM system.'''
        from numpy import vectorize
        from .constants import atomic_number
        atm_num = vectorize(atomic_number)
        return atm_num(self.atoms)

    @property
    def dim_atomic_numbers(self):
        '''Returns the atomic number of each atom in the QM system.'''
        from numpy import vectorize
        from .constants import atomic_number
        atm_num = vectorize(atomic_number)
        return atm_num(self.dim_atoms)


    @property
    def dim_masses(self):
        '''Returns the mass of each atom in DIM system.'''
        from numpy import vectorize
        from .constants import atomic_mass

        # Make the atomic mass function broadcastable on a numpy array
        atm_mass = vectorize(atomic_mass)
        # Find mass for each atom
        return atm_mass(self.dim_atoms)


    @property
    def allmasses(self):
        '''Returns the mass of each atom in entire system.'''
        from numpy import vectorize
        from .constants import atomic_mass

        # Make the atomic mass function broadcastable on a numpy array
        atm_mass = vectorize(atomic_mass)
        # Find mass for each atom
        return atm_mass(self.allatoms)


    @property
    def molecular_mass(self):
        '''Determines the molecular mass of all QM atoms.'''
        # Sum masses and return
        return self.masses.sum()


    @property
    def dim_mass(self):
        '''Determines the molecular mass of all DIM atoms.'''
        # Sum masses and return
        return self.dim_masses.sum()


    @property
    def system_mass(self):
        '''Determines the molecular mass of ALL atoms.'''
        # Sum masses and return
        return self.allmasses.sum()


    @property
    def center_of_mass(self):
        '''Determines the molecular center of mass of the QM atoms.'''
        from numpy import array
        # Transpose the coordinates for easy multiplicaiton
        C = self.coordinates.transpose()
        M = self.masses / self.masses.sum()
        return array([ (M*C[0]).sum(), (M*C[1]).sum(), (M*C[2]).sum() ])


    @property
    def center_of_nuc_charge(self):
        '''Determines the center of nuclear charge of the QM atoms.'''
        from numpy import array
        C = self.coordinates.transpose()
        M = self.atomic_numbers / self.atomic_numbers.sum()
        return array([ (M*C[0]).sum(), (M*C[1]).sum(), (M*C[2]).sum() ])


    @property
    def dim_center_of_mass(self):
        '''Determines the center of mass of the DIM atoms.'''
        from numpy import array
        # Only return if dim atoms and coordinates are not None
        if self.dim_atoms is not None and self.dim_coordinates is not None:
            C = self.dim_coordinates.transpose()
            M = self.dim_masses / self.dim_mass
            return array([ (M*C[0]).sum(), (M*C[1]).sum(), (M*C[2]).sum() ])
        else:
            return None


    def make_dim_atoms(self, a1=1, a2=None):
        '''Takes the specified atoms in the coordinates array and transfers
        them to the :py:attr:`dim_coordinates` array.  The atoms chosen are
        specified by *a1* and *a2*.

        If there are already :py:attr:`dim_atoms` present, the new atoms will
        be appended to the front.

        '''
        from numpy import delete

        # Convert atom number to index
        a1 -= 1
        # Default to all.
        if a2 is None: a2 = self.natoms

        # If DIM atoms are present, append chosen atoms to the front.
        # Otherwise, just copy the coordinates
        if 'DIM' in self.calctype:
            from numpy import concatenate
            self.dim_coordinates = concatentate((self.coordinates[a1:a2],
                                                 self.dim_coordinates))
            self.dim_atoms = concatentate((self.atoms[a1:a2], self.dim_atoms))
        else:
            self.dim_coordinates = self.coordinates[a1:a2].copy()
            self.dim_atoms = self.atoms[a1:a2].copy()
        self.ndim = len(self.dim_atoms)
        self.dim_elements = set(self.dim_atoms)
        self.ndim_elements = len(self.dim_elements)
        self.calctype.add('DIM')

        # Remove the atoms from the original coordinates
        self.coordinates = delete(self.coordinates, slice(a1, a2), axis=0)
        self.atoms = delete(self.atoms, slice(a1, a2))
        self.natoms = len(self.atoms)
        self.elements = set(self.atoms)
        self.nelements = len(self.elements)


    def find_center(self, type='geometrical', qm=True, dim=False):
        '''Locates the center of the molecule.

        There are two types of centers you can locate:

         - 'geometrical'
         - 'center-of-mass'

        The *dim* and *qm* booleans control what set of coordinates
        are examined.

        Raises :py:exc:`AssertionError` if *dim* is True and not a
        DIMQM calculation, if an invalid mode is entered, or if file
        is 'xyz' and the filename attribute does not exist

        '''
        if dim: assert 'DIM' in self.calctype, ('find_center(): '
                                                    'Not a DIMQM calculation.')

        if type == 'center-of-mass':
            wt = self.masses
        elif type == 'geometrical':
            wt = None
        else:
            raise ChemDataError ('find_center(): Invalid type: '+type)

        # Now find the center
        from numpy import average
        if dim and qm:
            return average(self.allcoords, axis=0, weights=wt)
        elif dim:
            return average(self.dim_coordinates, axis=0, weights=wt)
        else:
            return average(self.coordinates, axis=0, weights=wt)


    def shift_to_origin(self, type='geometrical'):
        '''Shifts the coordinates so that the center of the molecule
        is at the origin.

        *type* determines the center type.  See :py:func:`find_center`.

        This will shift ALL coordinates.  If this is a DIM calculation,
        the DIM system will be taken into account in the shifting.

        '''
        try:
            center = self.find_center(type, dim=True)
        except AssertionError:
            center = self.find_center(type)
        self.translate_coordinates(-center)


    def radii(self, set='vis'):
        '''Returns the radii for each :py:attr:`atom <atoms>` in the
        QM system.

        The option *set* specifies from which set of data you wish
        to collect the radii.  The options are:

            -vis: These radii are good for molecular visualizations but are not physical
            -vdw: These use the van Der Waals radii and are intended to be physical, however not all elements are available.

        '''
        from numpy import zeros
        from .constants import atomic_radius

        # Find radii for each atom
        rad = zeros(self.natoms)
        for i in range(self.natoms):
            rad[i] = atomic_radius(self.atoms[i], set)
        return rad


    def dim_radii(self, set='vis'):
        '''Same as :py:attr:`radii` but for :py:attr:`DIM atoms <dim_atoms>`'''
        from numpy import zeros
        from .constants import atomic_radius
        if 'COORDDEPEND' in self.key:
            return self.dim_cnradii
        rad = zeros(self.ndim)
        for i in range(self.ndim):
            rad[i] = atomic_radius(self.dim_atoms[i], set)
        return rad


    def allradii(self, set='vis'):
        '''Same as :py:attr:`radii` but for :py:attr:`all atoms <allatoms>`'''
        from numpy import zeros
        from .constants import atomic_radius

        rad = zeros(self.nallatoms)
        for i in range(self.nallatoms):
            rad[i] = atomic_radius(self.allatoms[i], set)
        return rad


    @property
    def bonds(self):
        '''Given that N is the number of bonds, this function returns
        an N x 2 numpy array containing the index of atom 1 and atom 2 of
        each bond.  Only calculates bonds for QM atoms.

        Therefore, bond i goes from self.coordinates[self.bonds[i,0]]
        to self.coordinates[self.bonds[i,1]].

        '''
        from numpy import array, asarray
        from .f2py import calc_bonds

        # Use info from input if alfready determind for us
        if self.program == 'ADF' and 'GUIBONDS' in self.key:
            bonds = []
            for g in self.key['GUIBONDS']:
                bonds.append(g.split()[1:3])
            # Subtract 1 to make indices then return
            return array(bonds, dtype=int) - 1

        # Otherwise figure it out ourselves.
        else:
            # Use a FORTRAN subroutine to find bonds.
            # Because FORTRAN uses fixed-length arrays, we return the number of
            # bonds found and use this to truncate the array
            try:
                bm, nbond = calc_bonds(self.coordinates.T, self.radii(), 1.1)
            except AttributeError:
                return None
            else:
                return asarray(bm.T[0:nbond], dtype=int)


    @property
    def dim_bonds(self):
        '''Given that N is the number of bonds, this function returns
        an N x 2 numpy array containing the index of atom 1 and atom 2 of
        each bond.  Only calculates bonds for DIM atoms

        Therefore, bond i goes from self.dim_coordinates[self.bonds[i,0]]
        to self.dim_coordinates[self.bonds[i,1]].

        '''
        from numpy import asarray
        from .f2py import calc_bonds
        assert 'DIM' in self.calctype, ('dim_bonds(): '
                                                  'Not a DIM calculation.')
        # Use a FORTRAN subroutine to find bonds.
        # Because FORTRAN uses fixed-length arrays, we return the number of
        # bonds found and use this to truncate the array
        try:
            bm, nbond = calc_bonds(self.dim_coordinates.T, self.dim_radii[0], 1.1)
        except AttributeError:
            return None
        else:
            return asarray(bm.T[0:nbond], dtype=int)


    def maxdist(self, dim=False, qm=True):
        '''Return the maximum distance between atoms in the system.

        By default, only the :py:attr:`QM atoms <atoms>` are inspected, but this can be
        edited with the *dim* and *qm*  booleans.

        '''
        from .f2py import minmax_pdist
        if dim: assert 'DIM' in self.calctype, ('maxdist(): '
                                                  'Not a DIM calculation.')
        if dim and qm:
            return minmax_pdist(self.allcoords.T)[1]
        elif dim:
            return minmax_pdist(self.dim_coordinates.T)[1]
        else:
            return minmax_pdist(self.coordinates.T)[1]


    def mindist(self, dimqm=False):
        '''Return the mimimum distance between atoms in the system.

        The boolean 'dimqm' determines if the min distance will be a
        distance between a DIM atom and a QM atom, or if it will be
        between any two atoms.

        Raises :py:exc:`AssertionError` if *dimqm* is True an this is not
        a DIMQM calculation.

        '''
        if dimqm: assert 'DIM' in self.calctype, ('mindist(): '
                                                  'Not a DIM calculation.')
        if dimqm:
            # Grab the subset where atom i is QM and atom j is DIM
            from .f2py import minmax_cdist
            return minmax_cdist(self.coordinates.T, self.dim_coordinates.T)[0]
        else:
            from .f2py import minmax_pdist
            return minmax_pdist(self.allcoords.T)[0]


    def translate_coordinates(self, transvec):
        '''Translates the coordinates given some translation vector,
        *transvec*, which is a numpy 1x3 array.

        '''
        from copy import copy
        from .f2py import translate
        tempvec = copy(transvec)
        if self.coordinates is not None:
            self.coordinates = translate(tempvec, self.coordinates)
        if 'DIM' in self.calctype:
            self.dim_coordinates = translate(tempvec, self.dim_coordinates)

    def rotate_coordinates(self, rotmat=None, angle=None, dir=None, rad=False,
        about='origin', properties=True):
        '''Rotates the coordinates given some rotation matrix, *rotmat*
        or angle, *angle*.

        May be rotated about the origin, the center-of-mass (com) or the
        center-of-nuclear-charge (conc).

        If using a predefined rotation matrix, it must be a 3x3 numpy array.

        If using an angles, the X, Y, or Z direction must also be specified
        with *dir*, which is not case-sensitive.

        By default, the angle input unit is degrees.  This can be changed
        by setting the boolean *rad* to True.

        The *properties* keyword signals to rotate the molecular properties
        (such as normal modes, polarizability, hyperpolarizability) by the
        same rotation. This is true by default.

        '''
        from numpy import array, dot, eye, einsum, zeros

        assert rotmat is not None or angle is not None, (
                "rotate_coordinates(): Must choose one of 'rotmat' or 'angle'")
        assert not (rotmat is not None and angle is not None), (
           "rotate_coordinates(): Must choose only one of 'rotmat' or 'angle'")

        # Create matrix if not given explicitly
        # Create a general 3-D rotation
        if rotmat is None:
            from math import radians, cos, sin
            assert dir is not None, (
                  "rotate_coordinates(): 'dir' must not be empty with 'angle'")
            assert ((isinstance(angle, list) and isinstance(dir, list)) or
                not (isinstance(angle, list) and isinstance(dir, list))), (
                "rotate_coordinates(): 'angle' and 'dir' must both be either "
                "a 'list' or 'str'")
            # Make sure lists are the same length, and turn non-list in to list
            if isinstance(angle, list):
                assert len(angle) == len(dir), ("rotate_coordinates(): "
                                   "'angle' and 'dir' must be the same length")
            else:
                angle = [angle]
                dir = [dir]

            # Generate the general rotation matrix
            rotmat = eye(3)
            for a, d in zip(angle, dir):
                if not rad: a = radians(a)
                c = cos(a)
                s = sin(a)
                if d.lower() == 'x':
                    temp = array([[ 1,  0,  0 ],
                                  [ 0,  c, -s ],
                                  [ 0,  s,  c ]], dtype=float)
                elif d.lower() == 'y':
                    temp = array([[ c,  0,  s ],
                                  [ 0,  1,  0 ],
                                  [-s,  0,  c ]], dtype=float)
                elif d.lower() == 'z':
                    temp = array([[ c, -s,  0 ],
                                  [ s,  c,  0 ],
                                  [ 0,  0,  1 ]], dtype=float)
                else:
                    raise ValueError ("Unknown value for 'dir': "+str(d))

                # Create the rotmat in rotation order. Not that this is
                # backwards from what you would expect because internally the
                # rotation is done with the coordinates first.
                rotmat = dot(temp, rotmat)

        # Define desired origin.
        if about.lower() == 'com':
            origin = self.center_of_mass
        elif about.lower() == 'conc':
            origin = self.center_of_nuc_charge
        else:
            origin = zeros((3))

        # Rotate. (NB: numpy.einsum is sufficiently fast.)
        self.coordinates = self.coordinates - origin
        self.coordinates = einsum('ab,ib->ia', rotmat, self.coordinates)
        self.coordinates = self.coordinates + origin

        if 'DIM' in self.calctype:
            self.dim_coordinates = self.dim_coordinates - origin
            self.dim_coordinates = einsum('ab,ib->ia', rotmat, self.dim_coordinates)
            self.dim_coordinates = self.dim_coordinates + origin

        # Rotate polarizabilities and normal modes if requested
        if properties:
            self.rotate_polarizabilities(rotmat)
            if 'FREQUENCIES' in self.calctype:
                self.normal_modes = einsum('ab,ijb->ija', rotmat, self.normal_modes)
            if 'EXCITATIONS' in self.calctype:
                self.TDM = einsum('ab,ib->ia', rotmat, self.TDM)


    def order_coords(self, atom=None, coord=None, dim=False, qm=True):
        '''Reorders the coordinates according to proximity to either
        a specific atom or a point in space.

        The *dim* and *qm* booleans control what set of coordinates
        are ordered.

        *atom* is the atom to reorder by.  If both *qm* and *dimqm* are
        True, the counting starts with the QM atoms.

        *coord* is the point in space to reorder with respect to as a.
        1 x 3 numpy array or a list.

        Both *atom* and *coord* cannot be true and raises
        :py:exc:`AssertionError` if they are.

        Raises :py:exc:`AssertionError` if *dim* is True and this is not a
        DIMQM calculation.

        '''
        from numpy import argsort, array
        from .f2py import calc_dist

        assert not (atom is None and coord is None), ('order_coords():'
                                        ' Must choose one of atoms or coords.')
        assert (atom is not None or coord is not None), ('order_coords(): '
                                           'Cannot use both atoms and coords.')
        if dim: assert 'DIM' in self.calctype, ('order_coords(): '
                                                       'Not a DIM calculation')

        # First determine the distance from the points
        if atom is not None:
            if dim and qm:
                if atom < self.natoms:
                    coord = self.coordinates[atom-1]
                else:
                    coord = self.dim_coordinates[self.natoms+atom-1]
            elif dim:
                coord = self.dim_coordinates[atom-1]
            else:
                coord = self.coordinates[atom-1]

        if qm:
            # Use Fortran routine calc_dist to determine distances from point.
            # Sort the distances and return the sorted indices.
            index = argsort(calc_dist(self.coordinates.T, coord))
            # Sort atoms coodinates and modes.
            self.atoms = self.atoms[index]
            self.coordinates = self.coordinates[index]
            if 'FREQUENCIES' in self.calctype:
                self.normal_modes = self.normal_modes[index]

        # Repeat for DIM
        if dim:
            index = argsort(calc_dist(self.dim_coordinates.T, coord))
            self.dim_atoms = self.dim_atoms[index]
            self.dim_coordinates = self.dim_coordinates[index]


    def step_size(self, sR=0.01):
        '''Return the mass weighted step-size for each normal mode.

        The option *sR*, which is the unweighted step size that was
        used to make the plus and minus direction coordinates for the
        3-point numerical differentiation, is defaulted to 0.01.

        If no normal modes are collected, :py:exc:`AssertionError` is raised.

        '''
        from numpy import zeros_like, array, sqrt, hstack
        from numpy.linalg import norm

        # The notation and equations are taken from Z. Phys. Chem., 217, 91.
        # See also Seth Morton's thesis, chapter 3, section 3.2.
        # R = mode coordinates
        # Q = mass-weighted mode coordinates
        # Qnorm = normalized mass-weighted coordinates
        # Qmag = |Q|
        # mass = atomic masses
        # rt_mass = square root of the atomic masses
        # Rmag = |R|
        # sR = coordinate stepsize
        # sQ = mass-weighted stepsize
        #
        # The stepsize is 0.01/|R|, but instead of calculating |R| directly,
        # they use Q and divide out the weights.  We will also do this to
        # maintain clarity, even though it adds in the extra step of dividing
        # out the weights.

        assert 'FREQUENCIES' in self.calctype, ('stepsize(): '
                                         'No normal modes have been collected')

        Q = self.normal_modes.copy() # Copy so we don't clobber original
        Rmag = zeros_like(self.v_frequencies)

        # Make the square root of atomic mass array
        # the same shape as the normal modes
        mass = array([self.masses]).T
        rt_mass = hstack((sqrt(mass), sqrt(mass), sqrt(mass)))

        for i in range(self.nmodes):
            Q[i] *= rt_mass
            Qnorm = Q[i] / norm(Q[i])
            # Find root of sum of squares to get the magnitude
            Rmag[i] = norm(Qnorm / rt_mass)

        return sR / Rmag

    def modes_to_tcl(self, freqmin=400., freqmax=2000., scalefactor=1.,
        vectorwidth=3., scalevector=1.0, color = "blue"):
        '''Creates .tcl files with normal mode vectors that may be easily
        viewed using VMD.'''

        assert 'FREQUENCIES' in self.calctype, ('modes_to_tcl(): '
                             'No normal modes have been collected')

        # Allow for degenerate modes
        deg_list = ('', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')
        degeneracy = 0
        previous_mode = ''

        # Cycle over each normal mode
        for im in range(self.nmodes):

            freq = self.v_frequencies[im] * scalefactor

            # Make sure that this mode is within range
            if freq < freqmin: continue
            if freq > freqmax: continue

            # Check for degeneracy
            strmode = '{0:.2f}'.format(freq)
            if previous_mode == strmode:
                degeneracy += 1
                strmode = strmode+'_'+deg_list[degeneracy]
            else:
                degeneracy = 0

            # Open filename for this mode
            fw = open('mode{0}-vmd.tcl'.format(strmode), 'w')
            previous_mode = strmode

            # Cycle over each atom
            print('draw color {0}'.format(color),file=fw)
            for ia in range(self.natoms):

                print('vmd_draw_vector2 0 {0} {2: 11.7f} {3: 11.7f} {4: 11.7f}'
                  '{1} {0} {5: 11.7f} {6: 11.7f} {7: 11.7f} {1} {8:4.2f}'.format(
                  '{', '}', self.coordinates[ia][0], self.coordinates[ia][1],
                  self.coordinates[ia][2], self.normal_modes[im][ia][0]*scalevector,
                  self.normal_modes[im][ia][1]*scalevector, self.normal_modes[im][ia][2]*scalevector,
                  vectorwidth), file=fw)

            # Close file
            fw.close()

    def modes_to_pymol(self, freqmin=400., freqmax=2000., scalefactor=1.,
        vectorwidth=0.2, scalevector=1.0, color = "blue",component='all'):
        '''Creates .pymol files with normal mode vectors that may be easily
        viewed using pymol.

        requires running the cgo_arrow.py user script before loading
        these files in pymol'''

        assert 'FREQUENCIES' in self.calctype, ('modes_to_pymol(): '
                             'No normal modes have been collected')

        # Allow for degenerate modes
        deg_list = ('', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i')
        degeneracy = 0
        previous_mode = ''

        # Cycle over each normal mode
        for im in range(self.nmodes):

            freq = self.v_frequencies[im] * scalefactor

            # Make sure that this mode is within range
            if freq < freqmin: continue
            if freq > freqmax: continue

            # Check for degeneracy
            strmode = '{0:.2f}'.format(freq)
            if previous_mode == strmode:
                degeneracy += 1
                strmode = strmode+'_'+deg_list[degeneracy]
            else:
                degeneracy = 0

            # Open filename for this mode
            fw = open('mode{0}.pymol'.format(strmode), 'w')
            previous_mode = strmode

            # Cycle over each atom
            for ia in range(self.natoms):

                if component == 'all':
                    print('cgo_modevec( {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f}, color="{9}")'.format(
                  '[', ']', self.coordinates[ia][0], self.coordinates[ia][1],
                  self.coordinates[ia][2], self.normal_modes[im][ia][0]*scalevector,
                  self.normal_modes[im][ia][1]*scalevector, self.normal_modes[im][ia][2]*scalevector,
                  vectorwidth, color), file=fw)
                if component == 'x':
                    print('cgo_modevec( {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f})'.format(
                  '[', ']', self.coordinates[ia][0], self.coordinates[ia][1],
                  self.coordinates[ia][2], self.normal_modes[im][ia][0]*scalevector,
                  0.0000, 0.0000,
                  vectorwidth), file=fw)
                if component == 'y':
                    print('cgo_modevec( {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f})'.format(
                  '[', ']', self.coordinates[ia][0], self.coordinates[ia][1],
                  self.coordinates[ia][2], 0.00000,
                  self.normal_modes[im][ia][1]*scalevector, 0.00000,
                  vectorwidth), file=fw)
                if component == 'z':
                    print('cgo_modevec( {0}{2: 11.7f},{3: 11.7f},{4: 11.7f}'
                  '{1}, {0}{5: 11.7f},{6: 11.7f},{7: 11.7f}{1}, radius={8:4.2f})'.format(
                  '[', ']', self.coordinates[ia][0], self.coordinates[ia][1],
                  self.coordinates[ia][2], 0.00000,
                  0.00000, self.normal_modes[im][ia][2]*scalevector,
                  vectorwidth), file=fw)

            # Close file
            fw.close()

    def atomic_dipoles_to_tcl(self, scalefactor=0.1, vectorwidth=3.):
        '''Creates aa atomic_dipoles.tcl file with atomic dipole vectors,
        which may be opened using VMD.'''

        assert self.atomic_multipole_moments is not None, ('No atomic multipole moments collected.')

        vec = self.atomic_multipole_moments['dipole'][:] * scalefactor
        fw = open('atomic_dipoles.tcl', 'w')

        for ia in range(self.natoms):
            print('vmd_draw_vector2 0 {0} {2: 11.7f} {3: 11.7f} {4: 11.7f}'
                  '{1} {0} {5: 11.7f} {6: 11.7f} {7: 11.7f} {1} {8:4.2f}'.format(
                  '{', '}', self.coordinates[ia][0], self.coordinates[ia][1],
                  self.coordinates[ia][2], vec[ia][0], vec[ia][1], vec[ia][2],
                  vectorwidth), file=fw)

        fw.close()

    def dtoLine(self, pointA, pointB, pointC):
        '''
        calculate the distance of pointC from an infinite line
        between pointB and pointA

        all pointA, pointB, and pointC are numpy arrays of shape 3
        returns distance

        $latex \frac{| \vec{AC} \times \vec{AB} |}{| \vec{AB} |} $
        this gives the height of a parallelogram created by sides AB and AC
        this height, is a perpendicular line to the original line that passes through AB
        '''
        import numpy as np

        vecAC = pointC - pointA
        vecAB = pointB - pointA

        num = np.cross(vecAC, vecAB)
        denum = vecAB

        # need magnitude of our vectors
        #use square root of dot product of itself
        num = np.sqrt(num.dot(num))
        denum = np.sqrt(denum.dot(denum))

        answer = num / denum
        distance = abs(answer)

        return distance

    def cutCylinder(self, atomA, atomB, radius, dim=False):
        '''
        atomA and atomB are in python counting scheme starting at 0
        calculates cylinder of the total system that has the dimensions
        Length = atomB - atomA
        radius = radius
        and finds all atoms that exist within that cylinder.
        Returns arrays of index and coordinates of the atoms in the cylinder

        NOTE: it currently calculates assuming the points A and B are perfectly aligned
                on the Z z-axis. Should be generalized for any system later
        '''
        import numpy as np

        index = []
        cylinder = []
        atomType = []
        #get points that will be used to create line
        if not dim:
            pointi = self.coordinates[atomA,:]
            pointj = self.coordinates[atomB,:]
            coords = self.coordinates
            atoms = self.atoms
        else:
            pointi = self.dim_coordinates[atomA,:]
            pointj = self.dim_coordinates[atomB,:]
            coords = self.dim_coordinates
            atoms = self.dim_atoms
        
        for i, value in enumerate(atoms):
            dist = self.dtoLine(pointi, pointj, coords[i,:])
            if dist <= radius:
                if pointi[2] <= coords[i,2] <= pointj[2]:
                    cylinder.append([coords[i,0], coords[i,1], coords[i,2]])
                    index.append(i)
                    atomType.append(atoms[i])
        cylinder = np.asarray(cylinder)
        index = np.asarray(index)
        atomType = np.asarray(atomType)
        return index, atomType, cylinder


    def cutWatCylinder(self, atomA, atomB, radius, dim=False):
        '''
        atomA and atomB are in python counting scheme starting at 0
        calculates cylinder of the total system that has the dimensions
        Length = atomB - atomA
        radius = radius
        and finds all atoms that exist within that cylinder.
        Returns arrays of index and coordinates of the atoms in the cylinder

        NOTE: it currently calculates assuming the points A and B are perfectly aligned
                on the Z z-axis. Should be generalized for any system later
                This version of the function works specifically for the solvated
                nanoparticle junction where DIM system has xyz coordinates in order
                O H H O H H ... for all 45000 water atoms that come first
        '''
        import numpy as np

        index = []
        cylinder = []
        atomType = []
        #get points that will be used to create line
        if not dim:
            pointi = self.coordinates[atomA,:]
            pointj = self.coordinates[atomB,:]
            coords = self.coordinates
            atoms = self.atoms
        else:
            pointi = self.dim_coordinates[atomA,:]
            pointj = self.dim_coordinates[atomB,:]
            coords = self.dim_coordinates
            atoms = self.dim_atoms
        for i, value in enumerate(atoms):
            if i % 3 != 0:
                continue
            if i > 45000:
                continue
            dist1 = self.dtoLine(pointi, pointj, coords[i,:])
            dist2 = self.dtoLine(pointi, pointj, coords[i+1,:])
            dist3 = self.dtoLine(pointi, pointj, coords[i+2,:])
            if (dist1 or dist2 or dist3) <= radius:
                if pointi[2] <= coords[i,2] <= pointj[2]:
                    cylinder.append([coords[i,0], coords[i,1], coords[i,2]])
                    cylinder.append([coords[i+1,0], coords[i+1,1], coords[i+1,2]])
                    cylinder.append([coords[i+2,0], coords[i+2,1], coords[i+2,2]])
                    index.append(i)
                    index.append(i+1)
                    index.append(i+2)
                    atomType.append(atoms[i])
                    atomType.append(atoms[i+1])
                    atomType.append(atoms[i+2])
        cylinder = np.asarray(cylinder)
        index = np.asarray(index)
        atomType = np.asarray(atomType)
        return index, atomType, cylinder
