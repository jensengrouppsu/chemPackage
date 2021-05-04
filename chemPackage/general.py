from __future__ import print_function, division
from .errorclass import ChemDataError

class General(object):
    '''Extends the ChemData class with methods to manipulate general
    properites.
    '''

    def printDipole(self, mag=False, dim=False, qm=True):
        '''Pretty print the dipole moment vector to standard output.

        :param dim: Print the :py:attr:`DIM dipole moment <chem.dim_dipole_tot>` vector.
        :param qm: Print the :py:attr:`QM dipole moment <chem.dipole>` vector.
        :param mag: Print the magnitude of the dipole moment vector.

        If both *qm* and *dim* are True, the dipole moments of
        the two systems are summed together.

        '''
        from numpy import vdot, sqrt

        assert self.dipole is not None, ('printDipole(): '
                                         'No dipole moment was collected.')
        if dim: assert 'DIM' in self.calctype, ('printDipole(): '
                                                    'Not a DIMQM calculation.')

        if dim and qm:
            dip = self.dipole + self.dim_dipole_tot
        elif dim:
            dip = self.dim_dipole_tot
        else:
            dip = self.dipole

        print()
        print('{0:>11}{1:>11}{2:>11}'.format('X', 'Y', 'Z'))
        print('{0[0]:14.4f}{0[1]:11.4f}{0[2]:11.4f}'.format(dip))
        if mag:
            print()
            print('MAGNITUDE: {0:11.6f}'.format(sqrt(vdot(dip, dip))))
        print()

    def printDipder(self):

        assert self.dgdip is not None, ('printDipole(): '
                                         'No dipole derivatives were collected.')

        print()
        print('{0:>11}{1:>11}{2:>11}{3:>11}'.format('Mode', 'X', 'Y', 'Z'))
        for x in range(self.nmodes):
            tmp = self.v_frequencies[x]
            print('{0:13.4f}'.format(tmp),'{0[0]:11.4f}{0[1]:11.4f}{0[2]:11.4f}'.format(self.dgdip[x]))

    def printPolarizability(self):
        '''Pretty-print the polarizability tensor to standard output.'''

        if self.program == 'DIM':
            assert self.dim_pol is not None, ('printPolarizability(): '
                                                 'No polarizability collected')
            a = self.dim_pol
        else:
            assert self.polarizability is not None, ('printPolarizability(): '
                                                 'No polarizability collected')
            a = self.polarizability
        X = ['X','Y','Z']
        for k in range(len(a)):
            imag = isinstance(a[k][0][0],complex)
            print()
            print('Polarizability Tensor: ')
            if imag: print ('Real:')
            print('{0:>9}{1:>11}{2:>11}'.format('X','Y','Z'))
            for i in range(3):
                print('{0:>1}{1:11.4f}{2:11.4f}{3:11.4f}'.format(X[i],
                      a[k][i][0].real, a[k][i][1].real, a[k][i][2].real))
            if imag:
                print ('Imaginary:')
                print('{0:>9}{1:>11}{2:>11}'.format('X','Y','Z'))
                for i in range(3):
                    print('{0:>1}{1:11.4f}{2:11.4f}{3:11.4f}'.format(X[i],
                          a[k][i][0].imag, a[k][i][1].imag, a[k][i][2].imag))
            print()

    def printSecondhyperpolarizability(self):
        '''Pretty-print the negative part of the second hyperpolarizability tensor to standard output.'''

        assert self.secondhyperpolarizability is not None,('printSecondHyperpolarizability(): '
                                                            'No second hyperpolarizability collected')
        rshp = self.secondhyperpolarizability.real
        ishp = self.secondhyperpolarizability.imag
        for k in range(len(rshp)):
            print()
            print('The "Negative Part" of Second Hyperpolarizability Tensor: ')
            label = '{0:>39}{1:30}{2:>16}'
            head  = '{0:>11}{1:>12}{2:>12}{3:8}{0:>11}{1:>12}{2:>12}'
            fr    = '{0:>4}{1:12.4f}{2:12.4f}{3:12.4f}  '
            fi    = '{0:>4}{1:12.4f}{2:12.4f}{3:12.4f}'
            print(label.format('Real','' ,'Imaginary'))
            print('              ', head.format('X', 'Y', 'Z', ''))
            print('           ',
                  fr.format('X', rshp[k,0,0,0,0], rshp[k,0,0,0,1], rshp[k,0,0,0,2]),
                  fi.format('X', ishp[k,0,0,0,0], ishp[k,0,0,0,1], ishp[k,0,0,0,2]))
            print('         X ',
                  fr.format('Y', rshp[k,0,0,1,0], rshp[k,0,0,1,1], rshp[k,0,0,1,2]),
                  fi.format('Y', ishp[k,0,0,1,0], ishp[k,0,0,1,1], ishp[k,0,0,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,0,0,2,0], rshp[k,0,0,2,1], rshp[k,0,0,2,2]),
                  fi.format('Z', ishp[k,0,0,2,0], ishp[k,0,0,2,1], ishp[k,0,0,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,0,1,0,0], rshp[k,0,1,0,1], rshp[k,0,1,0,2]),
                  fi.format('X', ishp[k,0,1,0,0], ishp[k,0,1,0,1], ishp[k,0,1,0,2]))
            print('   X     Y ',
                  fr.format('Y', rshp[k,0,1,1,0], rshp[k,0,1,1,1], rshp[k,0,1,1,2]),
                  fi.format('Y', ishp[k,0,1,1,0], ishp[k,0,1,1,1], ishp[k,0,1,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,0,1,2,0], rshp[k,0,1,2,1], rshp[k,0,1,2,2]),
                  fi.format('Z', ishp[k,0,1,2,0], ishp[k,0,1,2,1], ishp[k,0,1,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,0,2,0,0], rshp[k,0,2,0,1], rshp[k,0,2,0,2]),
                  fi.format('X', ishp[k,0,2,0,0], ishp[k,0,2,0,1], ishp[k,0,2,0,2]))
            print('         Z ',
                  fr.format('Y', rshp[k,0,2,1,0], rshp[k,0,2,1,1], rshp[k,0,2,1,2]),
                  fi.format('Y', ishp[k,0,2,1,0], ishp[k,0,2,1,1], ishp[k,0,2,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,0,2,2,0], rshp[k,0,2,2,1], rshp[k,0,2,2,2]),
                  fi.format('Z', ishp[k,0,2,2,0], ishp[k,0,2,2,1], ishp[k,0,2,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,1,0,0,0], rshp[k,1,0,0,1], rshp[k,1,0,0,2]),
                  fi.format('X', ishp[k,1,0,0,0], ishp[k,1,0,0,1], ishp[k,1,0,0,2]))
            print('         X ',
                  fr.format('Y', rshp[k,1,0,1,0], rshp[k,1,0,1,1], rshp[k,1,0,1,2]),
                  fi.format('Y', ishp[k,1,0,1,0], ishp[k,1,0,1,1], ishp[k,1,0,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,1,0,2,0], rshp[k,1,0,2,1], rshp[k,1,0,2,2]),
                  fi.format('Z', ishp[k,1,0,2,0], ishp[k,1,0,2,1], ishp[k,1,0,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,1,1,0,0], rshp[k,1,1,0,1], rshp[k,1,1,0,2]),
                  fi.format('X', ishp[k,1,1,0,0], ishp[k,1,1,0,1], ishp[k,1,1,0,2]))
            print('   Y     Y ',
                  fr.format('Y', rshp[k,1,1,1,0], rshp[k,1,1,1,1], rshp[k,1,1,1,2]),
                  fi.format('Y', ishp[k,1,1,1,0], ishp[k,1,1,1,1], ishp[k,1,1,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,1,1,2,0], rshp[k,1,1,2,1], rshp[k,1,1,2,2]),
                  fi.format('Z', ishp[k,1,1,2,0], ishp[k,1,1,2,1], ishp[k,1,1,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,1,2,0,0], rshp[k,1,2,0,1], rshp[k,1,2,0,2]),
                  fi.format('X', ishp[k,1,2,0,0], ishp[k,1,2,0,1], ishp[k,1,2,0,2]))
            print('         Z ',
                  fr.format('Y', rshp[k,1,2,1,0], rshp[k,1,2,1,1], rshp[k,1,2,1,2]),
                  fi.format('Y', ishp[k,1,2,1,0], ishp[k,1,2,1,1], ishp[k,1,2,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,1,2,2,0], rshp[k,1,2,2,1], rshp[k,1,2,2,2]),
                  fi.format('Z', ishp[k,1,2,2,0], ishp[k,1,2,2,1], ishp[k,1,2,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,2,0,0,0], rshp[k,2,0,0,1], rshp[k,2,0,0,2]),
                  fi.format('X', ishp[k,2,0,0,0], ishp[k,2,0,0,1], ishp[k,2,0,0,2]))
            print('         X ',
                  fr.format('Y', rshp[k,2,0,1,0], rshp[k,2,0,1,1], rshp[k,2,0,1,2]),
                  fi.format('Y', ishp[k,2,0,1,0], ishp[k,2,0,1,1], ishp[k,2,0,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,2,0,2,0], rshp[k,2,0,2,1], rshp[k,2,0,2,2]),
                  fi.format('Z', ishp[k,2,0,2,0], ishp[k,2,0,2,1], ishp[k,2,0,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,2,1,0,0], rshp[k,2,1,0,1], rshp[k,2,1,0,2]),
                  fi.format('X', ishp[k,2,1,0,0], ishp[k,2,1,0,1], ishp[k,2,1,0,2]))
            print('   Z     Y ',
                  fr.format('Y', rshp[k,2,1,1,0], rshp[k,2,1,1,1], rshp[k,2,1,1,2]),
                  fi.format('Y', ishp[k,2,1,1,0], ishp[k,2,1,1,1], ishp[k,2,1,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,2,1,2,0], rshp[k,2,1,2,1], rshp[k,2,1,2,2]),
                  fi.format('Z', ishp[k,2,1,2,0], ishp[k,2,1,2,1], ishp[k,2,1,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[k,2,2,0,0], rshp[k,2,2,0,1], rshp[k,2,2,0,2]),
                  fi.format('X', ishp[k,2,2,0,0], ishp[k,2,2,0,1], ishp[k,2,2,0,2]))
            print('         Z ',
                  fr.format('Y', rshp[k,2,2,1,0], rshp[k,2,2,1,1], rshp[k,2,2,1,2]),
                  fi.format('Y', ishp[k,2,2,1,0], ishp[k,2,2,1,1], ishp[k,2,2,1,2]))
            print('           ',
                  fr.format('Z', rshp[k,2,2,2,0], rshp[k,2,2,2,1], rshp[k,2,2,2,2]),
                  fi.format('Z', ishp[k,2,2,2,0], ishp[k,2,2,2,1], ishp[k,2,2,2,2]))
            print()

    def printOptRot(self):
        '''Pretty-print the optical rotation tensor to standard output.'''

        assert self.ord is not None, ('printOptRot(): '
                                      'No optical rotation tensor collected')
        a = self.ord
        X = ['X','Y','Z']
        for k in range(len(a)):
            imag = isinstance(a[k][0][0],complex)
            print()
            print('Optical Rotation (-G/omega) Tensor: ')
            if imag: print ('Real:')
            print('{0:>9}{1:>11}{2:>11}'.format('X','Y','Z'))
            for i in range(3):
                print('{0:>1}{1:11.4f}{2:11.4f}{3:11.4f}'.format(X[i],
                      a[k][i][0].real, a[k][i][1].real, a[k][i][2].real))
            if imag:
                print ('Imaginary:')
                print('{0:>9}{1:>11}{2:>11}'.format('X','Y','Z'))
                for i in range(3):
                    print('{0:>1}{1:11.4f}{2:11.4f}{3:11.4f}'.format(X[i],
                          a[k][i][0].imag, a[k][i][1].imag, a[k][i][2].imag))
            print()

    def printGtensor(self, gstensor=False):
        '''Pretty-print the G-tensor to standard output.'''

        assert self.gtensor is not None, ('printGtensor(): '
                                          'No G-tensor collected')
        if gstensor:
            assert self.gstensor is not None, ('printGtensor(): '
                                          'No Gs-tensor collected')
        if gstensor:
            a = self.gstensor
        else:
            a = self.gtensor
        X = ['X','Y','Z']
        for k in range(len(a)):
            imag = isinstance(a[k][0][0],complex)
            print()
            if gstensor:
                print('Magnetic Dipole-Electric Dipole (G) Tensor: ')
            else:
                print('Electric Dipole-Magnetic Dipole (G) Tensor: ')
            if imag: print ('Real:')
            print('{0:>9}{1:>11}{2:>11}'.format('X','Y','Z'))
            for i in range(3):
                print('{0:>1}{1:11.4f}{2:11.4f}{3:11.4f}'.format(X[i],
                      a[k][i][0].real, a[k][i][1].real, a[k][i][2].real))
            if imag:
                print ('Imaginary:')
                print('{0:>9}{1:>11}{2:>11}'.format('X','Y','Z'))
                for i in range(3):
                    print('{0:>1}{1:11.4f}{2:11.4f}{3:11.4f}'.format(X[i],
                          a[k][i][0].imag, a[k][i][1].imag, a[k][i][2].imag))
            print()

    def printAtensor(self, astensor=False):
        '''Pretty-print the A-tensor to standard output.'''

        assert self.atensor is not None, ('printAtensor(): '
                                          'No A-tensor collected')
        if astensor:
            assert self.astensor is not None, ('printAtensor(): '
                                          'No As-tensor collected')
            self.fold_Atensor('FOLD', astensor=True)
        else:
            self.fold_Atensor('FOLD')
        for k in range(len(self.atensor)):
            if astensor:
                A = self.astensor[k][:]
            else:
                A = self.atensor[k][:]
            imag = False
            if isinstance(A[0][0], complex):
                imag = True
            B = ['X','Y','Z']
            print()
            if astensor:
                print('Electric Quadrupole-Dipole (A) Tensor:')
            else:
                print('Electric Dipole-Quadrupole (A) Tensor:')
            if imag: print ('Real:')
            # print real part
            print('{0:>15}{1:>16}{2:>16}{3:>16}{4:>16}{5:>16}'.format('XX','XY',
                  'XZ','YY','YZ','ZZ'))
            for i in range(3):
                print('{0:>1}{1:16.4f}{2:16.4f}{3:16.4f}{4:16.4f}{5:16.4f}'
                      '{6:16.4f}'.format(B[i], A[i][0].real, A[i][1].real,
                      A[i][2].real, A[i][3].real, A[i][4].real, A[i][5].real))
            if imag:
                # print imagiary part
                print ('Imaginary:')
                print('{0:>15}{1:>16}{2:>16}{3:>16}{4:>16}{5:>16}'.format('XX','XY',
                      'XZ','YY','YZ','ZZ'))
                for i in range(3):
                    print('{0:>1}{1:16.4f}{2:16.4f}{3:16.4f}{4:16.4f}{5:16.4f}'
                          '{6:16.4f}'.format(B[i], A[i][0].imag, A[i][1].imag,
                          A[i][2].imag, A[i][3].imag, A[i][4].imag, A[i][5].imag))
            print()
            # Consider the trace and print necessary warnings
            trace = [0, 0, 0]
            for i in range(3):
                trace[i] = A[i][0] + A[i][3] + A[i][5]
                if abs(trace[i]/A[i].sum()) > 1e-3:
                    print('Warning, '+B[i]+'-component not traceless: '+str(trace[i]))

    def printCtensor(self):
        '''Pretty-print the C-tensor to standard output.'''

        assert self.ctensor is not None, ('printCtensor():'
                                          'No C-tensor collected')
        self.fold_Ctensor('FOLD')
        for k in range(len(self.ctensor)):
            C = self.ctensor[k][:]
            imag = False
            if isinstance(C[0][0], complex): imag = True
            S = ['XX','XY','XZ','YY','YZ','ZZ']
            print ()
            print ('Electric Quadrupole-Quadrupole (C) Tensor:')
            if imag: print ('Real:')
            print('{0:>16}{1:>16}{2:>16}{3:>16}{4:>16}{5:>16}'.format('XX','XY',
                  'XZ','YY','YZ','ZZ'))
            for i in range(6):
                print('{0:>2}{1:16.4f}{2:16.4f}{3:16.4f}{4:16.4f}{5:16.4f}'
                      '{6:16.4f}'.format(S[i], C[i][0].real, C[i][1].real,
                      C[i][2].real, C[i][3].real, C[i][4].real, C[i][5].real))
            if imag:
                print ('Imaginary:')
                print('{0:>16}{1:>16}{2:>16}{3:>16}{4:>16}{5:>16}'.format('XX','XY',
                      'XZ','YY','YZ','ZZ'))
                for i in range(6):
                    print('{0:>2}{1:16.4f}{2:16.4f}{3:16.4f}{4:16.4f}{5:16.4f}'
                          '{6:16.4f}'.format(S[i], C[i][0].imag, C[i][1].imag,
                          C[i][2].imag, C[i][3].imag, C[i][4].imag, C[i][5].imag))
            print()
            # Consider the trace and print necessary warnings
            for i in range(6):
                trace  = abs(C[i][0] + C[i][3] + C[i][5])
                tmax = ( abs(C[i]) ).max()
                if abs(trace/max(C[i].max(), abs(C[i].min()))) > 1e-3:
                    print('Warning, '+S[i]+'-component not traceless: '+str(trace))

    def printDtensor(self, dstensor=False):
        '''Pretty-print the D-tensor to standard output.'''

        assert self.dtensor is not None, ('printDtensor():'
                                          'No D-tensor collected')
        if dstensor:
            assert self.dstensor is not None, ('printDtensor():'
                                          'No Ds-tensor collected')
        self.fold_Dtensor('FOLD')
        for k in range(len(self.dtensor)):
            if dstensor:
                D = self.dstensor[k][:]
            else:
                D = self.dtensor[k][:]
            imag = False
            if isinstance(D[0][0], complex): imag = True
            S = ['X', 'Y', 'Z']
            print ()
            if dstensor:
                print ('Magnetic Dipole-Electric Quadrupole (D) Tensor:')
                D = D.transpose()
            else:
                print ('Electric Quadrupole-Magnetic Dipole (D) Tensor:')
            if imag: print ('Real:')
            print('{0:>15}{1:>16}{2:>16}{3:>16}{4:>16}{5:>16}'.format('XX','XY',
                  'XZ','YY','YZ','ZZ'))
            for i in range(3):
                print('{0:>1}{1:16.4f}{2:16.4f}{3:16.4f}{4:16.4f}{5:16.4f}'
                      '{6:16.4f}'.format(S[i], D[0][i].real, D[1][i].real,
                      D[2][i].real, D[3][i].real, D[4][i].real, D[5][i].real))
            if imag:
                print ('Imaginary:')
                print('{0:>16}{1:>16}{2:>16}{3:>16}{4:>16}{5:>16}'.format('XX','XY',
                      'XZ','YY','YZ','ZZ'))
                for i in range(3):
                    print('{0:>2}{1:16.4f}{2:16.4f}{3:16.4f}{4:16.4f}{5:16.4f}'
                          '{6:16.4f}'.format(S[i], D[0][i].imag, D[1][i].imag,
                          D[2][i].imag, D[3][i].imag, D[4][i].imag, D[5][i].imag))
            print()
            # Consider the trace and print necessary warnings
            for i in range(3):
                trace  = D[0][i] + D[3][i] + D[5][i]
                if abs(trace) > D.max() / 1000.:
                    print('Warning, '+S[i]+'-component not traceless: '+str(trace))

    def printQuadrupole(self):
        '''Pretty-print the quadrupole moment to standard output.'''

        assert self.quadrupole is not None, ('printQuadrupole(): '
                                              'No quadrupole moment collected')

        B = ['XX','XY','XZ','YY','YZ','ZZ']
        print ()
        print ('Quadrupole Moment (Buckingham Convention):')
        for i in range(6):
            print (B[i]+': {0:11.4f}'.format(self.quadrupole[i]))
        print ()


    def printEnergy(self):
        '''Pretty-print the various energy terms of the system to
        standard output.

        .. seealso::
         :py:attr:`~.pauli`
         | :py:attr:`~.electrostatic`
         | :py:attr:`~.steric`
         | :py:attr:`~.orbital`
         | :py:attr:`~.one-electron`
         | :py:attr:`~.two-electron`
         | :py:attr:`~.Coulomb`
         | :py:attr:`~.XC`
         | :py:attr:`~.nuc. repulsion`
         | :py:attr:`~.total`
         | :py:attr:`~.total HF`
         | :py:attr:`~.total DFT`
         | :py:attr:`~.MP2 correlation`
         | :py:attr:`~.CISD correlation`
         | :py:attr:`~.CCSD correlation`
         | :py:attr:`~.total MP2`
         | :py:attr:`~.total CISD`
         | :py:attr:`~.total CCSD`

        '''
        from .constants import HART2EV as H2EV

        assert isinstance(self.energy, dict), ('printEnergy(): '
                                        'No energies collected')

        head = {
                'pauli':            'Pauli Repulsion',
                'electrostatic':    'Electrostatic Interaction',
                'steric':           'Steric Interaction',
                'orbital':          'Orbital Interactions',
                'one-electron':     'Core Hamiltonian',
                'two-electron':     'Coulomb + Exchange',
                'Coulomb':          'Coulomb Interaction',
                'XC':               'Exchange-Correlation',
                'nuc. repulsion':   'Nuclear Repulsion',
                'total':            'Total',
                'total HF':         'Total HF',
                'total DFT':        'Total DFT',
                'MP2 correlation':  'MP2 Correlation',
                'CISD correlation': 'CISD Correlation',
                'CCSD correlation': 'CCSD Correlation',
                'total MP2':        'Total HF + MP2',
                'total CISD':       'Total HF + CISD',
                'total CCSD':       'Total HF + CCSD',
               }

        headlist = ('pauli', 'electrostatic', 'steric', 'orbital',
                    'one-electron', 'two-electron', 'Coulomb', 'XC',
                    'nuc. repulsion', 'total', 'total HF', 'total DFT',
                    'MP2 correlation', 'total MP2', 'CISD correlation',
                    'total CISD', 'CCSD correlation', 'total CCSD',)

        # Format string
        f = '{0:26}: {1:>12.3f} eV  => {2:>16.8f} a.u.'
        # Format string for orbital energies
        fo = '    {0:22}: {1:>9.3f} eV  => {2:>14.8f} a.u.'

        # Print each energy term
        print()
        for h in headlist:
            if h not in self.energy: continue
            if h == 'orbital':
                print('Orbital Interactions')
                for o in self.energy[h]:
                    if o == 'total': continue
                    print(fo.format(o, self.energy[h][o]*H2EV,
                                             self.energy[h][o]))
                print(fo.format('Total', self.energy[h]['total']*H2EV,
                                         self.energy[h]['total']))
            else:
                print(f.format(head[h], self.energy[h]*H2EV, self.energy[h]))
        print()

        # optional DIM term
        if 'DIM' in self.calctype:
            print(f.format('DIM/QM Interaction', self.dimqm_energy*H2EV,
                                                 self.dimqm_energy))
            if self.dim_energy is not None:
                print(f.format('Isolated DIM', self.dim_energy*H2EV,
                                               self.dim_energy))
            print()

    def cross_section(self, property=None, **kwargs):
        '''Returns the cross-section appropriate for this calculation type.

        This is just a wrapper that calls one of:
        :py:attr:`raman_cross_section`  (Raman)
        :py:attr:`absorption_cross_section` (Polarizability or Excitations)

        See the documentation of the called methods for the available options.

        More to be implemented soon!

        You may also specifically specify a cross_section type with
        the *property* keyword, i.e. for :py:attr:`raman_cross_section`
        *property* should be 'raman' and for
        :py:attr:`absorption_cross_section`  *property* should be 'absorption'.

        Raises an :py:exc:`AssertionError` if *property* is incorrect or if
        there is no calctype that fits.
        '''

        # Determine the method to call
        method = ''
        if property is None:
            if 'RAMAN' in self.calctype:
                method = 'raman'
            elif 'POLARIZABILITY' in self.calctype and 'FD' in self.calctype:
                method = 'absorption'
            elif 'EXCITATION' in self.calctype:
                method = 'absorption'
            else:
                raise AssertionError ('cross_section(): no calctype matches')
        elif property == 'raman':
            method = 'raman'
        elif property == 'absorption':
            method = 'absorption'
        else:
            raise AssertionError ('cross_section(): '
                                  'invalid property ({0})'.format(property))

        # Call the appropriate method
        if method == 'raman':
            return self.raman_cross_section(**kwargs)
        elif method == 'absorption':
            return self.absorption_cross_section(**kwargs)

    def absorb(self, other, ignore=set()):
        '''Method to "absorb" all of the data from another
        :py:class:`~.ChemData` instance to the current ChemData instance.

        :param other: The ChemData object to be absorbed.
        :type other: ChemObj
        :param ignore: Attributes that will **NOT** be absorbed.

        .. note:: Empty parameters in *other* will always be ignored

        '''
        from copy import deepcopy
        from .chemdata import ChemData
        # Not ChemData instance
        if not isinstance(other, ChemData):
            raise ChemDataError ('absorb(): other != ChemData()')
        # Make sure ignore is a set
        ignore = set(ignore)

        # Run over all attributes on the attribute list (minus the ones
        # in ignore) and copy into this instance.
        # Don't copy things that are None.
        for k in self._attrset.difference(ignore):
            if k not in other: continue
            setattr(self, k, deepcopy(getattr(other, k)))


    def empty(self, ignore=set()):
        '''Empties all attributes in the instance except those on the
        ignore list.

        '''

        # Make sure ignore is a set
        ignore = set(ignore)

        # Run over all attributes and empty them, except for the ones
        # in ignore.
        for k in self._attrset.difference(ignore):
            if k in ('calctype', 'subkey'):
                setattr(self, k, set())
            elif k == 'key':
                setattr(self, k, {})
            else:
                setattr(self, k, None)


    def include(self, other):
        '''Adds attributes from different calculations into the current
        instance.

        If both calculations have the same type of information (i.e.
        polarizabilities at the different frequencies) then the new data is
        appended, and sorted according to increasing energy.  If the data
        is the same (i.e. the same frequency) then it is not copied.

        Attributes added are:
          - :py:attr:`excitation_energies`
          - :py:attr:`oscillator_strengths`
          - :py:attr:`excitation_symmetries`
          - :py:attr:`excitation_type`
          - :py:attr:`TDM`
          - :py:attr:`transitions`
          - :py:attr:`v_frequencies`
          - :py:attr:`normal_modes`
          - :py:attr:`IR`
          - :py:attr:`_raman`
          - :py:attr:`deltas`
          - :py:attr:`e_frequencies`
          - :py:attr:`polarizability`
          - :py:attr:`dim_pol`
          - :py:attr:`dgdip`
          - :py:attr:`dtdip`

        All other attributes are simply ignored.

        '''
        from numpy import argsort, append, where
        from copy import copy, deepcopy

        # Copy excitations
        if 'EXCITATIONS' in other.calctype:

            # If both are same type, add unique values
            if 'EXCITATIONS' in self.calctype:

                # Convert transitions to list so we can append to it
                self.transitions = list(self.transitions)

                # Loop over all excitations to be added
                for i in range(other.nmodes):
                    # Next value if it already exists
                    if other.excitation_energies[i] in self.excitation_energies:
                        continue
                    # Otherwise, append the values
                    self.excitation_energies = append(
                                                 self.excitation_energies,
                                                [other.excitation_energies[i]],
                                                 axis=0)
                    self.oscillator_strengths = append(
                                                self.oscillator_strengths,
                                               [other.oscillator_strengths[i]],
                                                axis=0)
                    self.excitation_symmetries = append(
                                               self.excitation_symmetries,
                                              [other.excitation_symmetries[i]],
                                               axis=0)
                    self.excitation_types = append(self.excitation_types,
                                                  [other.excitation_types[i]],
                                                   axis=0)
                    self.TDM = append(self.TDM, [other.TDM[i]], axis=0)
                    self.transitions.append(other.transitions[i])
                    self.nexcite += 1

                # Sort the new values
                indx = argsort(self.v_frequencies)
                self.excitation_energies = self.excitation_energies[indx]
                self.oscillator_strengths = self.oscillator_strengths[indx]
                self.excitation_symmetries = self.excitation_symmetries[indx]
                self.excitation_type = self.excitation_type[indx]
                self.TDM = self.TDM[indx]
                temp = deepcopy(self.transitions)
                for i, ix in enumerate(indx):
                    self.transitions[i] = temp[ix]
                self.transitions = tuple(self.transitions)

            # Otherwise, simply copy information
            else:
                self.excitation_energies = copy(other.excitation_energies)
                self.oscillator_strengths = copy(other.oscillator_strengths)
                self.excitation_symmetries = copy(other.excitation_symmetries)
                self.excitation_types = copy(other.excitation_types)
                self.TDM = copy(other.TDM)
                self.dgdip = copy(other.dgdip)
                self.transitions = deepcopy(other.transitions)
                self.nexcite = len(self.excitation_energies)
                self.calctype.add('EXCITATIONS')

        # Copy normal modes
        if 'FREQUENCIES' in other.calctype:

            # If both are same type, add unique values
            if 'FREQUENCIES' in self.calctype:

                # Loop over all modes to be added
                for i in range(other.nmodes):
                    # Next value if it already exists
                    if other.v_frequencies[i] in self.v_frequencies:
                        continue
                    # Otherwise, append the values
                    self.normal_modes = append(self.normal_modes,
                                            [other.normal_modes[i]], axis=0)
                    self.v_frequencies = append(self.v_frequencies,
                                              [other.v_frequencies[i]], axis=0)
                    self.IR = append(self.IR, [other.IR[i]], axis=0)
                    if '_raman' in self and '_raman' in other:
                        self._raman = append(self._raman,
                                            [other._raman[i]], axis=0)
                    elif '_raman' in other:
                        self._raman = other._raman[i]
                    if 'deltas' in self and 'deltas' in other:
                        self.deltas = append(self.deltas,
                                            [other.deltas[i]], axis=0)
                    elif 'deltas' in other:
                        self.deltas = other.deltas[i]
                    self.nmodes += 1
                    if 'dgdip' in other:
                        self.dgdip = other.dgdip[i]
                    if 'dtdip' in other:
                        self.dtdip = other.dtdip[i]

                # Sort the new values
                indx = argsort(self.v_frequencies)
                self.v_frequencies = self.v_frequencies[indx]
                self.normal_modes = self.normal_modes[indx]
                self.IR = self.IR[indx]
                if '_raman' in self: self._raman = self._raman[indx]
                if 'deltas' in self: self.deltas = self.deltas[indx]

            # Otherwise, simply copy information
            else:
                self.v_frequencies = copy(other.v_frequencies)
                self.normal_modes = copy(other.normal_modes)
                self.IR = copy(other.IR)
                self._raman = copy(other._raman)
                self.deltas = copy(other.deltas)
                self.dgdip = copy(other.dgdip)
                self.dtdip = copy(other.dtdip)
                self.nmodes = len(self.v_frequencies)
                self.calctype.add('FREQUENCIES')
                if 'RAMAN' in other.calctype: self.calctype.add('RAMAN')
                if 'DELTAS' in other.calctype: self.calctype.add('DELTAS')

        # Copy polarizabilities
        if 'POLARIZABILITY' in other.calctype:

            # If both are same type, add unique values
            if 'POLARIZABILITY' in self.calctype:

                # Loop over all frequencies to be added
                for i in range(other.npol):
                    # Next value if it already exists
                    if other.e_frequencies[i] in self.e_frequencies:
                        continue
                    # Otherwise, append the values
                    try:
                        self.qm_pol = append(self.qm_pol,
                                            [other.qm_pol[i]], axis=0)
                    except TypeError:
                        pass
                    self.e_frequencies = append(self.e_frequencies,
                                              [other.e_frequencies[i]], axis=0)
                    if 'dim_pol' in self:
                        self.dim_pol = append(self.dim_pol,
                                                    [other.dim_pol[i]], axis=0)
                    if "dim_dipoles" in self:
                        self.dim_dipoles['FD scattered']=append(self.dim_dipoles['FD scattered'],[other.dim_dipoles['FD scattered'][i]],axis=0)
                    self.npol += 1
                    if 'dim_efficiencies' in self:
                        self.dim_efficiencies = append(self.dim_efficiencies, other.dim_efficiencies, axis=0)
                # Sort the new values
                indx = argsort(self.e_frequencies)
                self.e_frequencies = self.e_frequencies[indx]
                try:
                    self.qm_pol = self.qm_pol[indx]
                except TypeError:
                    pass
                if 'dim_pol' in self: self.dim_pol = self.dim_pol[indx]
                try:
                    self.dim_efficiencies = self.dim_efficiencies[indx]
                except TypeError:
                    passpr

            # Otherwise, simply copy information
            else:
                self.e_frequencies = copy(other.e_frequencies)
                self.qm_pol = copy(other.qm_pol)
                self.dim_pol = copy(other.dim_pol)
                self.npol = len(self.e_frequencies)
                self.atensor = copy(other.atensor)
                self.astensor = copy(other.astensor)
                self.gtensor = copy(other.gtensor)
                self.gstensor = copy(other.gstensor)
                self.dtensor = copy(other.dtensor)
                self.dstensor = copy(other.dstensor)
                self.ctensor = copy(other.ctensor)
                self.calctype.add('POLARIZABILITY')


    def collect_dir(self, dir=None, raise_err=True, files=[]):
        '''Collects from files in a directory and aggregates the results.

        The directory is assumed to be the current directory if none is given.
        If no filenames are given, the glob `*.out` is used to collect the
        files. Otherwise, only the listed files are collected, and globs are
        accepted.

        '''
        from . import collect
        from prep import natsort_key
        from numpy import append, array
        from glob import glob
        import os, sys

        # Copy all info from one file if this is an empty instance
        absorb_bool = True if not self else False

        # Assign correct directory
        if dir is None: dir = os.curdir
        # If files are given, assume each name is a glob and expand it.
        if files:
            tmp = []
            for f in files:
                tmp.extend(glob(os.path.join(dir, f)))
            files = [f for f in tmp]
        # Otherwise glob the directory
        else:
            files = glob(os.path.join(dir, '*.out'))

        # Sort the filenames 'naturally', i.e. in numerical order, not ASCII.
        files.sort(key=natsort_key)

        # Loop over the numbers and collect the info into a ChemData object
        self.e_frequencies = array([], dtype=float)
        for file in files:

            # Collect the data from this file
            try:
                d = collect(file, raise_err=raise_err)
            except IOError as e:
                sys.exit(str(e))

            # Include the information from the this file into self
            self.include(d)

            # For last file, absorb remaining information into self
            if file == files[-1] and absorb_bool:
                ignore = ('polarizability', 'qm_pol', 'dim_pol', 'e_frequencies', 'npol',
                          'excitation_energies', 'oscillator_strengths',
                          'excitation_symmetries', 'excitation_type', 'TDM',
                          'transitions', 'v_frequencies', 'normal_modes', 'IR',
                          '_raman', 'nmodes', 'nexcite', 'deltas')
                self.absorb(d, ignore)
                self.filename = os.path.join(dir, '*.out')
