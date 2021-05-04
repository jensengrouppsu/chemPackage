from __future__ import print_function, division
from .errorclass import ChemDataError

class Polarizability(object):
    '''Extends ChemData class with methods to manipulate polarizability.'''

    @property
    def polarizability(self):
        '''Smartly returns the polarizability of the system.

        If the calculation was a QM calculation, it will return :py:attr:`qm_pol`
        If the calculation was a DIM calculation, it will return :py:attr:`dim_pol`
        If the calculation was a DIM/QM calculation, it will return :py:attr:`qm_pol` + :py:attr:`dim_pol`
        '''
        try:
            p = self.qm_pol #+ self.dim_pol
            #print("QM + DIM")
        except TypeError:
            if self.qm_pol is not None:
                p = self.qm_pol
                #print("QM")
            else:
                p = self.dim_pol
                #print("DIM")
        return p

    def printTensor(self, property=None, iso=False, ani=False,
                          dim=False, qm=False, unit='au',
                          p1=1, p2=None):
        '''Pretty print the polarizability tensor to standard output.

        If the :py:attr:`polarizability tensor <polarizability>` is complex,
        then the real and imaginary tensors are printed alongside each other.

        The option *property* can be either 'pol' or 'ord' for polarizability
        or optical rotation, respectively.  If neither is given, the property
        one will be guessed based on the calctype.

        By default, only the tensor is printed.  You can add verbosity
        by requesting the :py:attr:`isotropic` or
        :py:attr:`anisotropic <anisotropic2>` polarizability to
        be printed with the *iso* and *ani* booleans.  Also, a header will
        be printed describing the tensor depending on the type of
        calculation.  If the calculation is **RAMAN** then the
        :py:attr:`scattering factor <scattering_factor>` is also printed.

        If the boolean *dim* is true, then the
        :py:attr:`DIM polarizability <dim_pol>` tensor will be printed.
        Same for the boolean *qm*. If both are true,
        the polarizabilities of the two systems are summed together.

        If this is not a RAMAN calculation, the
        :py:attr:`frequency <e_frequencies>` will be printed
        along with the tensor.  The argument *unit* allows you to choose
        the unit this is printed in,  either 'au' (Hartrees), 'eV',
        or 'nm'.

        *p1* is the lowest polarizability to print.  *p2* is the highest.
        This option allows you to print a subset of tensors collected.

        Raises :py:exc:`AssertionError` if no polarizabilities were collected.
        '''
        from numpy import array, sum
        import re

        if property == 'pol':
            assert 'POLARIZABILITY' in self.calctype, ('printTensor(): '
                                         'No polarizabilities were collected.')
        elif property == 'ord':
            assert 'OPTICAL ROTATION' in self.calctype, ('printTensor(): '
                                        'No optical rotations were collected.')
        elif property is None:
            # Try to guess the property
            if 'POLARIZABILITY' in self.calctype:
                property = 'pol'
            elif 'OPTICAL ROTATION' in self.calctype:
                property = 'ord'
            else:
                raise AssertionError ('printTensor(): Nothing to average.')

        if dim: assert 'DIM' in self.calctype, ('printTensor(): '
                                                    'Not a DIMQM calculation.')

        # Convert pol number to index
        p1 -= 1
        # Default to all.
        if p2 is None: p2 = self.npol

        # Set the unit and energy type
        if not re.match(r'au|ev|nm', unit, re.I):
            raise ChemDataError ('printTensor() : Invalid unit: '+unit)
        elif 'RAMAN' not in self.calctype:
            unit = { 'au' : 'a.u.', 'ev' : 'eV', 'nm' : 'nm' }[unit.lower()]
            etyp = { 'a.u.' : 'Frequency',
                     'eV'   : 'Energy',
                     'nm'   : 'Wavelength' }[unit]
        else:
            unit = 'cm-1'
            etyp = 'Normal Mode'

        # Set the tensor
        if property == 'pol':
            prop = 'Polarizability'
            if dim and qm:
                tensor = self.qm_pol + self.dim_pol
            elif dim:
                tensor = self.dim_pol
            # Added functionality to print qm_pol only. -Pengchong, Oct 2017
            elif qm:
                tensor = self.qm_pol 
            else:
                tensor = self.polarizability
        elif property == 'ord':
            prop = 'Optical Rotation'
            tensor = self.ord

        # Start the formatting of the label
        if 'RAMAN' in self.calctype:
           lbl = self.v_frequencies
        else:
           from .constants import HART2EV, HART2NM
           lbl = { 'a.u.' : self.e_frequencies,
                   'eV'   : HART2EV(self.e_frequencies),
                   'nm'   : HART2NM(self.e_frequencies)}[unit]

        # Diagonalized or not
        d = ', Diagonalized' if self._diagonalized else ''

        # Print selected tensors in class
        print()
        for n in range(p1, p2):

            t = tensor[n]
            r = t.real
            i = t.imag

            # Prints out the proper heading depending on the energy type
            # and the unit.
            # First, make the number at most 7 digits with 'good' formatting
            # Strip off whitespace
            s = '{0:7g}'.format(lbl[n]).strip()
            # Now place this number in the heading
            head = 'FD {4}: {0} {1} {2}{3}'.format(etyp, s, unit, d, prop)
            # Replace the heading if this was a static calculation
            if s == '0': head = 'Static {0}'.format(prop)
            print(head)

            # Pretty-print the tensor
            if 'FD' in self.calctype:
                label = '{0:>35}{1:35}{2:>22}'
                head  = '{0:>17}{1:>16}{2:>16}{3:6}{0:>17}{1:>16}{2:>16}'
                fr    = '{0:>4}{1[0]:16.3f}{1[1]:16.3f}{1[2]:16.3f}  '
                fi    = '{0:>4}{1[0]:16.3f}{1[1]:16.3f}{1[2]:16.3f}'
                print(label.format('Real', '', 'Imaginary'))
                print(head.format('X', 'Y', 'Z', ''))
                print(fr.format('X', r[0,:]), fi.format('X', i[0,:]))
                print(fr.format('Y', r[1,:]), fi.format('Y', i[1,:]))
                print(fr.format('Z', r[2,:]), fi.format('Z', i[2,:]))
            else:
                head = '{0:>17}{1:>17}{2:>17}'
                f    = '{0:>4}{1[0]:17.4f}{1[1]:17.4f}{1[2]:17.4f}'
                print(head.format('X', 'Y', 'Z'))
                print(f.format('X', r[0,:]))
                print(f.format('Y', r[1,:]))
                print(f.format('Z', r[2,:]))

            # Add the invariants at the end if requested
            if iso or ani or 'RAMAN' in self.calctype:
                print()
            if iso:
                print('ISOTROPIC PART    = ', end='')
                iso = self.isotropic(property=property, dim=dim, qm=qm)[n]
                if 'FD' in self.calctype:
                    print('{0:12g} + {1:7g}j'.format(iso.real, iso.imag))
                else:
                    print('{0:12g}'.format(iso))
            if ani:
                from math import sqrt
                ani = self.anisotropic2(property=property, dim=dim, qm=qm)[n]
                print('ANISOTROPIC PART  = {0:12g}'.format(sqrt(ani)))
            if 'RAMAN' in self.calctype:
                scat = self.scattering_factor()[n]
                print('SCATTERING FACTOR = {0:12g}'.format(scat))

            print()


    def printHyperpolarizability(self, unit='au'):
        '''Pretty-prints the
        :py:attr:`hyperpolarizability tensor <hyperpolarizability>`.

        Note that field one is the subtensor, field two is the row, and
        field three is the column.

        The argument *unit* allows you to choose the unit this is printed
        in,  either 'au' (Hartrees), 'eV', or 'nm'.

        Raises :py:exc:`AssertionError` if no hyperpolarizabilities were
        collected.
        '''
        import re

        assert 'HYPERPOLARIZABILITY' in self.calctype, (
             'printHyperpolarizability(): No hyperpolarizabilities collected.')

        if 'FD' in self.calctype:
            # Property type
            prop = 'Hyperpolarizability'

            hp = self.hyperpolarizability['FD']
            rhp = hp.real
            ihp = hp.imag

            #This "if statement" is set for ADF_FD. If there is no 'HYPERPOL' or
            #'TWONPLUSONE', jump to "head = " 
            if 'HYPERPOL' in self.subkey or 'TWONPLUSONE' in self.subkey:
                # Set the unit and energy type
                if not re.match(r'au|ev|nm', unit, re.I):
                    raise ChemDataError ('printTensor() : Invalid unit: '+unit)
                else:
                    unit = { 'au' : 'a.u.', 'ev' : 'eV', 'nm' : 'nm' }[unit.lower()]
                    etyp_b = { 'a.u.' : 'B-Freq.',
                               'eV'   : 'B-Energy',
                               'nm'   : 'B-Wavelength' }[unit]
                    etyp_c = { 'a.u.' : 'C-Freq.',
                               'eV'   : 'C-Energy',
                               'nm'   : 'C-Wavelength' }[unit]

                # Start the formatting of the label
                from .constants import HART2EV, HART2NM
                lbl_b = { 'a.u.' : self.b_e_frequencies,
                          'eV'   : HART2EV(self.b_e_frequencies),
                          'nm'   : HART2NM(self.b_e_frequencies)}[unit]
                lbl_c = { 'a.u.' : self.c_e_frequencies,
                          'eV'   : HART2EV(self.c_e_frequencies),
                          'nm'   : HART2NM(self.c_e_frequencies)}[unit]

                # Prints out the proper heading depending on the energy type 
                # and the unit.
                # First, make the number at most 7 digits with 'good' formatting
                # Strip off whitespace
                s_b = '{0:7g}'.format(lbl_b[0]).strip()
                s_c = '{0:7g}'.format(lbl_c[0]).strip()
                # Now place this number in the heading
                head = 'FD {6}: {0} {1} {2} {3} {4} {5}'.format(etyp_b, s_b, unit,
                                                                etyp_c, s_c, unit, prop)
                # Replace the heading if this was a static calculation
                if s_b == '0' and s_c == '0': head = 'Static {0}'.format(prop)
            else:
                # Finite Difference method can only do EOPE or STATIC calculation 
                head = 'EOPE Hyperpolarizability'
            print()
            print(head)

            # Pretty print the hyperpolarizability
            label = '{0:>33}{1:30}{2:>16}'
            head  = '{0:>11}{1:>12}{2:>12}{3:8}{0:>11}{1:>12}{2:>12}'
            fr    = '{0:>4}{1:12.3f}{2:12.3f}{3:12.3f}  '
            fi    = '{0:>4}{1:12.3f}{2:12.3f}{3:12.3f}'
            print(label.format('Real', '', 'Imaginary'))
            print('       ', head.format('X', 'Y', 'Z', ''))
            print('     ',
                  fr.format('X', rhp[0,0,0], rhp[0,0,1], rhp[0,0,2]),
                  fi.format('X', ihp[0,0,0], ihp[0,0,1], ihp[0,0,2]))
            print('   X ',
                  fr.format('Y', rhp[0,1,0], rhp[0,1,1], rhp[0,1,2]),
                  fi.format('Y', ihp[0,1,0], ihp[0,1,1], ihp[0,1,2]))
            print('     ',
                  fr.format('Z', rhp[0,2,0], rhp[0,2,1], rhp[0,2,2]),
                  fi.format('Z', ihp[0,2,0], ihp[0,2,1], ihp[0,2,2]))
            print()
            print('     ',
                  fr.format('X', rhp[1,0,0], rhp[1,0,1], rhp[1,0,2]),
                  fi.format('X', ihp[1,0,0], ihp[1,0,1], ihp[1,0,2]))
            print('   Y ',
                  fr.format('Y', rhp[1,1,0], rhp[1,1,1], rhp[1,1,2]),
                  fi.format('Y', ihp[1,1,0], ihp[1,1,1], ihp[1,1,2]))
            print('     ',
                  fr.format('Z', rhp[1,2,0], rhp[1,2,1], rhp[1,2,2]),
                  fi.format('Z', ihp[1,2,0], ihp[1,2,1], ihp[1,2,2]))
            print()
            print('     ',
                  fr.format('X', rhp[2,0,0], rhp[2,0,1], rhp[2,0,2]),
                  fi.format('X', ihp[2,0,0], ihp[2,0,1], ihp[2,0,2]))
            print('   Z ',
                  fr.format('Y', rhp[2,1,0], rhp[2,1,1], rhp[2,1,2]),
                  fi.format('Y', ihp[2,1,0], ihp[2,1,1], ihp[2,1,2]))
            print('     ',
                  fr.format('Z', rhp[2,2,0], rhp[2,2,1], rhp[2,2,2]),
                  fi.format('Z', ihp[2,2,0], ihp[2,2,1], ihp[2,2,2]))
            print()
        else:
            keys = ('STATIC', 'SHG', 'EOPE', 'OR')
            headers = tuple(['{0} Hyperpolarizability'.format(x) for x in keys])
            for key, head in zip(keys, headers):
                try:
# jbb5516 added .real to end of next line to take out imaginary parts
                    hp = self.hyperpolarizability[key].real
                except KeyError:
                    continue
                hp = hp.real
                print()
                print(head)
                print('     {0:>12}{1:>10}{2:>10}'.format('X', 'Y', 'Z'))
                print('     {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'X', hp[0,0,0], hp[0,0,1], hp[0,0,2]))
                print('  X  {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'Y', hp[0,1,0], hp[0,1,1], hp[0,1,2]))
                print('     {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'Z', hp[0,2,0], hp[0,2,1], hp[0,2,2]))
                print()
                print('     {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'X', hp[1,0,0], hp[1,0,1], hp[1,0,2]))
                print('  Y  {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'Y', hp[1,1,0], hp[1,1,1], hp[1,1,2]))
                print('     {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'Z', hp[1,2,0], hp[1,2,1], hp[1,2,2]))
                print()
                print('     {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'X', hp[2,0,0], hp[2,0,1], hp[2,0,2]))
                print('  Z  {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'Y', hp[2,1,0], hp[2,1,1], hp[2,1,2]))
                print('     {0:>4}{1:10.3f}{2:10.3f}{3:10.3f}'.format(
                            'Z', hp[2,2,0], hp[2,2,1], hp[2,2,2]))
                print()

    def printSecondHyperpolarizability(self, unit='au'):
        '''Pretty-prints the 
        :py:attr:`second hyperpolarizability tensor <second hyperpolarizability>`.

        Note that field one is the tensor, field two is the subtensor, field three
        is the row, and field four is the column.

        The argument *unit* allows you to choose the unit this is printed 
        in,  either 'au' (Hartrees), 'eV', or 'nm'.

        Raises :py:exc:`AssertionError` if no second hyperpolarizabilities were 
        collected.
        '''
        import re

        assert 'SECOND HYPERPOLARIZABILITY' in self.calctype, (
             'printSecondHyperpolarizability(): No second hyperpolarizabilities collected.')

        if 'FD' in self.calctype:
            # Property type
            prop = 'Second Hyperpolarizability'

            shp = self.secondhyperpolarizability['FD']
            rshp = shp.real
            ishp = shp.imag

            #This "if statement" is set for ADF_FD. If there is no 'GAMMA', jump to "head = " 
            if 'GAMMA' in self.subkey:
                # Set the unit and energy type
                if not re.match(r'au|ev|nm', unit, re.I):
                   raise ChemDataError ('printTensor() : Invalid unit: '+unit)
                else:
                        unit = { 'au' : 'a.u.', 'ev' : 'eV', 'nm' : 'nm'}[unit.lower()]
                        etyp_b = { 'a.u.' : 'B-Freq.',
                                   'eV'   : 'B-Energy',
                                   'nm'   : 'B-Wavelength' }[unit]
                        etyp_c = { 'a.u.' : 'C-Freq.',
                                   'eV'   : 'C-Energy',
                                   'nm'   : 'C-Wavelength' }[unit]
                        etyp_d = { 'a.u.' : 'D-Freq.',
                                   'eV'   : 'D-Energy',
                                   'nm'   : 'D-Wavelength' }[unit]

                # Start the formatting of the label
                from .constants import HART2EV, HART2NM
                lbl_b = { 'a.u.' : self.b_e_frequencies,
                          'eV'   : HART2EV(self.b_e_frequencies),
                          'nm'   : HART2NM(self.b_e_frequencies)}[unit]
                lbl_c = { 'a.u.' : self.c_e_frequencies,
                          'eV'   : HART2EV(self.c_e_frequencies),
                          'nm'   : HART2NM(self.c_e_frequencies)}[unit]
                lbl_d = { 'a.u.' : self.d_e_frequencies,
                          'eV'   : HART2EV(self.d_e_frequencies),
                          'nm'   : HART2NM(self.d_e_frequencies)}[unit]

                # Prints out the proper heading depending on the energy type
                # and the unit.
                # First, make the number at most 7 digits with 'good' formatting
                # Strip off whitespace
                s_b = '{0:7g}'.format(lbl_b[0]).strip()
                s_c = '{0:7g}'.format(lbl_c[0]).strip()
                s_d = '{0:7g}'.format(lbl_d[0]).strip()
                # Now place this number in the heading
                head = 'FD {9}: {0} {1} {2} {3} {4} {5} {6} {7} {8}'.format(etyp_b, s_b, unit,
                                                                            etyp_c, s_c, unit,
                                                                            etyp_d, s_d, unit, prop)
                # Replace the heading if this was a static calculation
                if s_b == '0' and s_c == '0' and s_d == '0': head = 'Static {0}'.format(prop)
            else:
                # Finite Difference method can only do EOPE or STATIC calculation 
                head = 'EFISHG Second Hyperpolarizability'
            print()
            print(head)

            # Pretty print the hyperpolarizability 
            label = '{0:>39}{1:30}{2:>16}'
            head  = '{0:>11}{1:>12}{2:>12}{3:8}{0:>11}{1:>12}{2:>12}'
            fr    = '{0:>4}{1:12.3f}{2:12.3f}{3:12.3f}  '
            fi    = '{0:>4}{1:12.3f}{2:12.3f}{3:12.3f}'
            print(label.format('Real','' ,'Imaginary'))
            print('              ', head.format('X', 'Y', 'Z', ''))
            print('           ',
                  fr.format('X', rshp[0,0,0,0], rshp[0,0,0,1], rshp[0,0,0,2]),
                  fi.format('X', ishp[0,0,0,0], ishp[0,0,0,1], ishp[0,0,0,2]))
            print('         X ',
                  fr.format('Y', rshp[0,0,1,0], rshp[0,0,1,1], rshp[0,0,1,2]),
                  fi.format('Y', ishp[0,0,1,0], ishp[0,0,1,1], ishp[0,0,1,2]))
            print('           ',
                  fr.format('Z', rshp[0,0,2,0], rshp[0,0,2,1], rshp[0,0,2,2]),
                  fi.format('Z', ishp[0,0,2,0], ishp[0,0,2,1], ishp[0,0,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[0,1,0,0], rshp[0,1,0,1], rshp[0,1,0,2]),
                  fi.format('X', ishp[0,1,0,0], ishp[0,1,0,1], ishp[0,1,0,2]))
            print('   X     Y ',
                  fr.format('Y', rshp[0,1,1,0], rshp[0,1,1,1], rshp[0,1,1,2]),
                  fi.format('Y', ishp[0,1,1,0], ishp[0,1,1,1], ishp[0,1,1,2]))
            print('           ',
                  fr.format('Z', rshp[0,1,2,0], rshp[0,1,2,1], rshp[0,1,2,2]),
                  fi.format('Z', ishp[0,1,2,0], ishp[0,1,2,1], ishp[0,1,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[0,2,0,0], rshp[0,2,0,1], rshp[0,2,0,2]),
                  fi.format('X', ishp[0,2,0,0], ishp[0,2,0,1], ishp[0,2,0,2]))
            print('         Z ',
                  fr.format('Y', rshp[0,2,1,0], rshp[0,2,1,1], rshp[0,2,1,2]),
                  fi.format('Y', ishp[0,2,1,0], ishp[0,2,1,1], ishp[0,2,1,2]))
            print('           ',
                  fr.format('Z', rshp[0,2,2,0], rshp[0,2,2,1], rshp[0,2,2,2]),
                  fi.format('Z', ishp[0,2,2,0], ishp[0,2,2,1], ishp[0,2,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[1,0,0,0], rshp[1,0,0,1], rshp[1,0,0,2]),
                  fi.format('X', ishp[1,0,0,0], ishp[1,0,0,1], ishp[1,0,0,2]))
            print('         X ',
                  fr.format('Y', rshp[1,0,1,0], rshp[1,0,1,1], rshp[1,0,1,2]),
                  fi.format('Y', ishp[1,0,1,0], ishp[1,0,1,1], ishp[1,0,1,2]))
            print('           ',
                  fr.format('Z', rshp[1,0,2,0], rshp[1,0,2,1], rshp[1,0,2,2]),
                  fi.format('Z', ishp[1,0,2,0], ishp[1,0,2,1], ishp[1,0,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[1,1,0,0], rshp[1,1,0,1], rshp[1,1,0,2]),
                  fi.format('X', ishp[1,1,0,0], ishp[1,1,0,1], ishp[1,1,0,2]))
            print('   Y     Y ',
                  fr.format('Y', rshp[1,1,1,0], rshp[1,1,1,1], rshp[1,1,1,2]),
                  fi.format('Y', ishp[1,1,1,0], ishp[1,1,1,1], ishp[1,1,1,2]))
            print('           ',
                  fr.format('Z', rshp[1,1,2,0], rshp[1,1,2,1], rshp[1,1,2,2]),
                  fi.format('Z', ishp[1,1,2,0], ishp[1,1,2,1], ishp[1,1,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[1,2,0,0], rshp[1,2,0,1], rshp[1,2,0,2]),
                  fi.format('X', ishp[1,2,0,0], ishp[1,2,0,1], ishp[1,2,0,2]))
            print('         Z ',
                  fr.format('Y', rshp[1,2,1,0], rshp[1,2,1,1], rshp[1,2,1,2]),
                  fi.format('Y', ishp[1,2,1,0], ishp[1,2,1,1], ishp[1,2,1,2]))
            print('           ',
                  fr.format('Z', rshp[1,2,2,0], rshp[1,2,2,1], rshp[1,2,2,2]),
                  fi.format('Z', ishp[1,2,2,0], ishp[1,2,2,1], ishp[1,2,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[2,0,0,0], rshp[2,0,0,1], rshp[2,0,0,2]),
                  fi.format('X', ishp[2,0,0,0], ishp[2,0,0,1], ishp[2,0,0,2]))
            print('         X ',
                  fr.format('Y', rshp[2,0,1,0], rshp[2,0,1,1], rshp[2,0,1,2]),
                  fi.format('Y', ishp[2,0,1,0], ishp[2,0,1,1], ishp[2,0,1,2]))
            print('           ',
                  fr.format('Z', rshp[2,0,2,0], rshp[2,0,2,1], rshp[2,0,2,2]),
                  fi.format('Z', ishp[2,0,2,0], ishp[2,0,2,1], ishp[2,0,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[2,1,0,0], rshp[2,1,0,1], rshp[2,1,0,2]),
                  fi.format('X', ishp[2,1,0,0], ishp[2,1,0,1], ishp[2,1,0,2]))
            print('   Z     Y ',
                  fr.format('Y', rshp[2,1,1,0], rshp[2,1,1,1], rshp[2,1,1,2]),
                  fi.format('Y', ishp[2,1,1,0], ishp[2,1,1,1], ishp[2,1,1,2]))
            print('           ',
                  fr.format('Z', rshp[2,1,2,0], rshp[2,1,2,1], rshp[2,1,2,2]),
                  fi.format('Z', ishp[2,1,2,0], ishp[2,1,2,1], ishp[2,1,2,2]))
            print()
            print('           ',
                  fr.format('X', rshp[2,2,0,0], rshp[2,2,0,1], rshp[2,2,0,2]),
                  fi.format('X', ishp[2,2,0,0], ishp[2,2,0,1], ishp[2,2,0,2]))
            print('         Z ',
                  fr.format('Y', rshp[2,2,1,0], rshp[2,2,1,1], rshp[2,2,1,2]),
                  fi.format('Y', ishp[2,2,1,0], ishp[2,2,1,1], ishp[2,2,1,2]))
            print('           ',
                  fr.format('Z', rshp[2,2,2,0], rshp[2,2,2,1], rshp[2,2,2,2]),
                  fi.format('Z', ishp[2,2,2,0], ishp[2,2,2,1], ishp[2,2,2,2]))
            print()
        else:
            keys = ('STATIC', 'THG', 'EFISHG', 'EFIOR', 'IDRI', 'OKE')
            headers = tuple(['{0} Second Hyperpolarizability'.format(x) for x in keys])
            for key, head in zip(keys, headers):
                try:
                    shp = self.secondhyperpolarizability[key]
                except KeyError:
                    continue
                rshp = shp.real
                fr    = '{0:>4}{1:12.3f}{2:12.3f}{3:12.3f}  '
                print()
                print(head)
                print('     {0:>21}{1:>12}{2:>12}'.format('X', 'Y', 'Z'))
                print('           ',
                      fr.format('X', rshp[0,0,0,0], rshp[0,0,0,1], rshp[0,0,0,2]))
                print('         X ',
                      fr.format('Y', rshp[0,0,1,0], rshp[0,0,1,1], rshp[0,0,1,2]))
                print('           ',
                      fr.format('Z', rshp[0,0,2,0], rshp[0,0,2,1], rshp[0,0,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[0,1,0,0], rshp[0,1,0,1], rshp[0,1,0,2]))
                print('   X     Y ',
                      fr.format('Y', rshp[0,1,1,0], rshp[0,1,1,1], rshp[0,1,1,2]))
                print('           ',
                      fr.format('Z', rshp[0,1,2,0], rshp[0,1,2,1], rshp[0,1,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[0,2,0,0], rshp[0,2,0,1], rshp[0,2,0,2]))
                print('         Z ',
                      fr.format('Y', rshp[0,2,1,0], rshp[0,2,1,1], rshp[0,2,1,2]))
                print('           ',
                      fr.format('Z', rshp[0,2,2,0], rshp[0,2,2,1], rshp[0,2,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[1,0,0,0], rshp[1,0,0,1], rshp[1,0,0,2]))
                print('         X ',
                      fr.format('Y', rshp[1,0,1,0], rshp[1,0,1,1], rshp[1,0,1,2]))
                print('           ',
                      fr.format('Z', rshp[1,0,2,0], rshp[1,0,2,1], rshp[1,0,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[1,1,0,0], rshp[1,1,0,1], rshp[1,1,0,2]))
                print('   Y     Y ',
                      fr.format('Y', rshp[1,1,1,0], rshp[1,1,1,1], rshp[1,1,1,2]))
                print('           ',
                      fr.format('Z', rshp[1,1,2,0], rshp[1,1,2,1], rshp[1,1,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[1,2,0,0], rshp[1,2,0,1], rshp[1,2,0,2]))
                print('         Z ',
                      fr.format('Y', rshp[1,2,1,0], rshp[1,2,1,1], rshp[1,2,1,2]))
                print('           ',
                      fr.format('Z', rshp[1,2,2,0], rshp[1,2,2,1], rshp[1,2,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[2,0,0,0], rshp[2,0,0,1], rshp[2,0,0,2]))
                print('         X ',
                      fr.format('Y', rshp[2,0,1,0], rshp[2,0,1,1], rshp[2,0,1,2]))
                print('           ',
                      fr.format('Z', rshp[2,0,2,0], rshp[2,0,2,1], rshp[2,0,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[2,1,0,0], rshp[2,1,0,1], rshp[2,1,0,2]))
                print('   Z     Y ',
                      fr.format('Y', rshp[2,1,1,0], rshp[2,1,1,1], rshp[2,1,1,2]))
                print('           ',
                      fr.format('Z', rshp[2,1,2,0], rshp[2,1,2,1], rshp[2,1,2,2]))
                print()
                print('           ',
                      fr.format('X', rshp[2,2,0,0], rshp[2,2,0,1], rshp[2,2,0,2]))
                print('         Z ',
                      fr.format('Y', rshp[2,2,1,0], rshp[2,2,1,1], rshp[2,2,1,2]))
                print('           ',
                      fr.format('Z', rshp[2,2,2,0], rshp[2,2,2,1], rshp[2,2,2,2]))
                print()

    def printOptical(self, property='absorption', qm=True, dim=False,
                     pathlength=None, concentration=None, unit='au',
                     absunit='angstroms'):
        '''Prints a list of optical properties.

        The choices of properties are 'absorption', 'absorptivitty',
        'absorbance', and 'transmission'.  'absorbance' and 'transmission'
        are calculated from beer's law and thus the pathlength and
        concentration are required.  'absorption' may be returned in a variety
        of units; see the :py:attr:`absorption_cross_section` documentation for
        the available options and this is controlled here with the absunit
        keyword.  The unit keyword controlls the frequency unit as in the
        :py:attr:`printTensor` method.
        '''
        import re
        from .constants import HART2EV, HART2NM

        # Prep for printing
        if property == 'absorption':
            values = self.absorption_cross_section(qm=qm, dim=dim, unit=unit)
            label = 'Absorption Cross-Section'
            if absunit == 'angstroms':
                u = unicode(u'\u212B^2/molecule')
            elif absunit == 'bohr':
                u = 'bohr^2/molecule'
            elif absunit == 'nm':
                u = 'nm^2/molecule'
            elif absunit == 'cm':
                u = 'cm^2/molecule'
            elif absunit == 'm':
                u = 'm^2/molecule'
        elif property == 'absorptivitty':
            values = self.molar_absorptivitty(qm=qm, dim=dim)
            label, u = 'Molar Absorptivitty', 'L mol^{-1} cm^{-1}'
        elif property == 'absorbance':
            values = self.absorbance(pathlength, concentration, qm=qm, dim=dim)
            label, u  = 'Absorbance', 'unitless'
        elif property == 'transmittance':
            values = self.transmittance(pathlength, concentration, qm=qm, dim=dim)
            label, u = 'Transmittance', 'unitless'
        else:
            raise ChemDataError ('printOptical(): '
                                 'invalid property ('+property+')')

        # Choose frequency unit
        if not re.match(r'au|ev|nm', unit, re.I):
            raise ChemDataError ('printTensor() : Invalid unit: '+unit)
        unit = { 'au' : 'a.u.', 'ev' : 'eV', 'nm' : 'nm' }[unit.lower()]
        etyp = { 'a.u.' : 'Frequency',
                 'eV'   : 'Energy',
                 'nm'   : 'Wavelength' }[unit]
        frequencies = { 'a.u.' : self.e_frequencies,
                        'eV'   : HART2EV(self.e_frequencies),
                        'nm'   : HART2NM(self.e_frequencies)}[unit]


        # Title
        if property in ('absorbance', 'transmittance'):
            string = 'Pathlength={0:.2f} cm, Concentration={1:.3E} M'
            print(string.format(pathlength, concentration))
        string = unicode(u'{0:>10} {1:<6} {2:>24} {3:<20}')
        print(string.format(etyp, '('+unit+')', label, '('+u+')'))

        # A format string
        if unit == 'a.u.':
            fmt = '{0:^17.4E} {1:^45.5E}'
        else:
            fmt = '{0:^17.4f} {1:^45.5E}'

        # Print the properties
        for freq, val in zip(frequencies, values):
            print(fmt.format(freq, val))


    @staticmethod
    def tensor_isotropic(tn):
        from numpy import trace
        return trace(tn, axis1=1, axis2=2) / 3


    @staticmethod
    def tensor_anisotropic2(tn):
        from numpy import empty, absolute
        # Calculate anisotripoc polarizability squared using eq.14 from
        # J. Chem. Phys, 75, 5615:
        # \frac{3}{4}\left[\sum_{ij}\alpha_{ij}\alpha_{ij)^\ast
        #                + \sum_{ij}\alpha_{ij}\alpha_{ji)^\ast\right]
        # - \frac{1}{2}\sum_{ij}\alpha_{ii}\alpha_{jj}^\ast

        # It should be noted that there is no complex anisotropic
        # polarizability like there is isotropic polarizability.  This is due
        # to the fact that the anisotropic polarizaibiltiy is instrinsically
        # a magnatude due to the squared factor, and as such must be real.
        cj = tn.conjugate()
        r = range(3)
        t = empty(len(tn), dtype=float)
        for n in range(len(tn)):
            cross1 = sum([tn[n,i,j] * cj[n,i,j] for i in r for j in r])
            cross2 = sum([tn[n,i,j] * cj[n,j,i] for i in r for j in r])
            diag   = sum([tn[n,i,i] * cj[n,j,j] for i in r for j in r])
            t[n] = absolute(( 3 / 4 ) * ( cross1 + cross2 ) - ( 1 / 2 ) * diag)
        return t


    def isotropic(self, property=None, qm=False, dim=False):
        '''Calculates the isotropic polarizability for the requested property
        and the requested system.

        The option *property* can be either 'pol' or 'ord' for polarizability
        or optical rotation, respectively.  If neither is given, the property
        one will be guessed based on the calctype.

        If the boolean *qm* is true, then only the isotropic values for the
        QM system are returned.  If *dim* is true, then the isotropic values
        for the DIM system will be returned.  If both are true, then the
        sum of the two are returned.  If both are false, it will try to return
        for the QM system first and if that is empty, then try the DIM system.
        '''
        if property == 'pol':
            assert 'POLARIZABILITY' in self.calctype, ('isotropic(): '
                                         'No polarizabilities were collected.')
        elif property == 'ord':
            assert 'OPTICAL ROTATION' in self.calctype, ('isotropic(): '
                                        'No optical rotations were collected.')
        elif property is None:
            # Try to guess the property
            if 'POLARIZABILITY' in self.calctype:
                property = 'pol'
            elif 'OPTICAL ROTATION' in self.calctype:
                property = 'ord'
            else:
                raise AssertionError ('isotropic(): Nothing to average.')
        if dim:
            assert 'DIM' in self.calctype, ('isotropic(): '
                                                    'Not a DIMQM calculation.')

        # Determine what tensor to return and do it
        if property == 'pol':
            if qm and dim:
                tensor = self.qm_pol + self.dim_pol
            elif dim:
                tensor = self.dim_pol
            elif qm:
                tensor = self.qm_pol
            else:
                tensor = self.polarizability
        else:
            tensor = self.ord
        return self.tensor_isotropic(tensor)


    def anisotropic2(self, property=None, qm=False, dim=False):
        '''Calculates the anisotropic polarizability squared
        The API is the same as for :py:attr:`isotropic`..'''

        if property == 'pol':
            assert 'POLARIZABILITY' in self.calctype, ('anisotropic2(): '
                                         'No polarizabilities were collected.')
        elif property == 'ord':
            assert 'OPTICAL ROTATION' in self.calctype, ('anisotropic2(): '
                                        'No optical rotations were collected.')
        elif property is None:
            # Try to guess the property
            if 'POLARIZABILITY' in self.calctype:
                property = 'pol'
            elif 'OPTICAL ROTATION' in self.calctype:
                property = 'ord'
            else:
                raise AssertionError ('anisotropic2(): Nothing to average.')
        if dim:
            assert 'DIM' in self.calctype, ('anisotropic2(): '
                                                    'Not a DIMQM calculation.')

        # Determine what tensor to return and do it
        if property == 'pol':
            if qm and dim:
                tensor = self.qm_pol + self.dim_pol
            elif dim:
                tensor = self.dim_pol
            elif qm:
                tensor = self.qm_pol
            else:
                tensor = self.polarizability

        else:
            tensor = self.ord
        return self.tensor_anisotropic2(tensor)


    def absorption_cross_section(self, unit='angstroms', qm=True, dim=False):
        '''Returns the absorption cross-section of the system.

        If the :py:attr:`calctype` is 'POLARIZABILITY', this will be calcualted
        from the imaginary component of the isotropic polarizability.
        If the :py:attr:`calctype` is 'EXCIATATION', this will be calculated
        from the excitation.

        By default, the cross-section is returned in units of
        angstroms^2/molecule.  This can be changed to nm^2/molecule,
        bohr^2/moleulce, cm^2/molecule, or m^2/molecule with the arguments
        'nm', 'bohr', 'cm', or 'm' to the keyword *unit*, respectively.

        QM and DIM may be toggled with the *qm* and *dim* keywords.

        A :py:exc:`AssertionError` is raised if the calctype does not
        allow this property to be calculated.
        '''
        from .constants import PI, LIGHT_AU, NOCONV, BOHR2ANGSTROM, BOHR2NM
        from .constants import BOHR2CM, BOHR2M

        # Assign the correct unit converter
        if unit.lower() == 'angstroms':
            conv = BOHR2ANGSTROM
        elif unit.lower() == 'bohr':
            conv = NOCONV
        elif unit.lower() == 'nm':
            conv = BOHR2NM
        elif unit.lower() == 'cm':
            conv = BOHR2CM
        elif unit.lower() == 'm':
            conv = BOHR2M
        else:
            raise AssertionError ('absorption_cross_section(): '
                                  'Invalid unit ('+unit.lower()+')')

        if 'POLARIZABILITY' in self.calctype and 'FD' in self.calctype:
            iso = self.isotropic(property='pol', qm=qm, dim=dim).imag
            # Conversion is performed twice because the property is squared
            return conv(conv(( 4 * PI * self.e_frequencies / LIGHT_AU ) * iso))
        elif 'EXCITATION' in self.calctype:
            pass # To be implemented
        else:
            raise AssertionError ('absorption_cross_section(): '
                          'Not a FD-POLARIZABILITY or EXCITATION calculation')


    def molar_absorptivitty(self, qm=True, dim=False):
        '''Calculates the molar absorptivity of the system.

        The unit is returned in L mol^{-1} cm^{-1} or M^{-1} cm^{-1}
        (these are equivalent).

        QM and DIM may be toggled with the *qm* and *dim* keywords.

        A :py:exc:`AssertionError` is raised if the calctype does not
        allow this property to be calculated.
        '''
        from math import log as ln
        from .constants import ANGSTROM2CM, AVOGADRO

        # Conversion factor to get from angstroms^2 to L/cm.
        # L == dm^3
        # dm^2 == 1E18 angstroms^2
        # dm^2 == 0.1 dm^3/cm
        # 0.1 dm^3/cm = 1E18 angstroms^2
        # 1E-19 = dm^3/(cm*angstroms^2)
        CONVFACTOR = 1E-19

        # Get the absorption cross-section in angstroms^2/molecule
        acs = self.absorption_cross_section(qm=qm, dim=dim)
        # Use this to calculate the molar absorptivity
        # LN(10) accounts for Beer's law
        return AVOGADRO * CONVFACTOR * acs / ln(10)


    def absorbance(self, pathlength, concentration, qm=True, dim=False):
        '''Calculates the absorbance of the system using Beer's law.

        The user must supply the *pathlength* in centimeters and the
        *concentration* in molarity.

        QM and DIM may be toggled with the *qm* and *dim* keywords.

        A :py:exc:`AssertionError` is raised if the calctype does not
        allow this property to be calculated.
        '''
        eps = self.molar_absorptivitty(qm=qm, dim=dim)
        return eps * pathlength * concentration


    def transmittance(self, pathlength, concentration, qm=True, dim=False):
        '''Calculates the transmittance from the absorbance.

        The transmittance is given on a scale from 0 to 1.  The user must
        convert to percent if she so desires.

        The user must supply the *pathlength* in centimeters and the
        *concentration* in molarity.

        QM and DIM may be toggled with the *qm* and *dim* keywords.

        A :py:exc:`AssertionError` is raised if the calctype does not
        allow this property to be calculated.
        '''
        from numpy import power
        absorb = self.absorbance(pathlength, concentration, qm=qm, dim=dim)
        return power(10, -absorb)


    def raman_depol_ratio(self):
        '''Calculates the Raman depolarization ratio.'''

        from numpy import zeros, empty, array, linalg

        assert 'POLARIZABILITY' in self.calctype, ('raman_depol_ratio(): '
                                     'No polarizabilities were collected.')

        # Arrays storing the total polarizability averages.
        if 'FREQUENCIES' in self.calctype:
            a_term = zeros((self.nmodes,10))
            a_sum  = empty((self.nmodes,3))
            raman_ratio = empty(self.nmodes)
            tn = self.qm_pol
            cj = tn.conjugate()
        else:
            print ('No frequency file is given')

        for n in range(len(tn)):
            a_term[n][0] = tn[n][0][0] + tn[n][1][1] + tn[n][2][2]
            a_term[n][1] = tn[n][0][1] - tn[n][1][0]
            a_term[n][2] = tn[n][0][2] - tn[n][2][0]
            a_term[n][3] = tn[n][1][2] - tn[n][2][1]
            a_term[n][4] = tn[n][0][1] + tn[n][1][0]
            a_term[n][5] = tn[n][0][2] + tn[n][2][0]
            a_term[n][6] = tn[n][1][2] + tn[n][2][1]
            a_term[n][7] = tn[n][0][0] - tn[n][1][1]
            a_term[n][8] = tn[n][0][0] - tn[n][2][2]
            a_term[n][9] = tn[n][1][1] - tn[n][2][2]

            a_sum[n][0] = 1/3 * linalg.norm(a_term[n][0])
            a_sum[n][1] = 1/2 * (linalg.norm(a_term[n][1]) + linalg.norm(a_term[n][2]) + linalg.norm(a_term[n][3])) 
            a_sum[n][2] = 1/2 * (linalg.norm(a_term[n][4]) + linalg.norm(a_term[n][5]) + linalg.norm(a_term[n][6]))
            + 1/3 * (linalg.norm(a_term[n][7]) + linalg.norm(a_term[n][8]) + linalg.norm(a_term[n][9]))

            raman_ratio[n] = (5*a_sum[n][1] + 3*a_sum[n][2]) / (10*a_sum[n][0] + 4*a_sum[n][2])

        return raman_ratio

    def hpol_average(self, ratio=False, real=False, imag=False, cv_mode=False, mode_bterm=False, component='all'):
        '''
        Calculates the average of the hyperpolarizability for hyper-Raman or hyper-Rayleigh
        scattering (HRS or HRayS) collected perpendicular to the incident light direction.
        This also allows for the depolarization ratio calculations for HRS at each mode or
        for HRayS at each frequency. Additionally, the contributions to HRays from both real
        and imaginary parts can be separately calculated and obtained, and for crystal violet,
        more specific information about its HRS can be obtained (hard coded at this moment).

        Additional options are:
       
        ratio      ==> Calcualte the depolarization ratio for either HRS or HRayS
        Real       ==> Calcluate the HRS by only using the real part of hyperpolarizability
        imag       ==> Calcluate the HRS by only using the imaginary part of hyperpolarizability
        cv_mode    ==> Calcluate all bterm contibutions (no coefficients involved) for the mode
                       at 1591.896 wavenumber^{-1} of CV
        mode_bterm ==> Calcluate components of the dominant bterm contributions (coeffcients
                       involved) to both <\\beta_ZZZ> and <\\beta_XZZ> for the mode at 1591.896
                       wavenumber^{-1} of CV
        '''

        from numpy import zeros, empty, array

        assert 'HYPERPOLARIZABILITY' in self.calctype, ('hpol_average(): '
                                'No hyperpolarizabilities were collected.')

        # Arrays storing the total hyperpolarizability averages.
        if 'FREQUENCIES' in self.calctype:
            b_avg = zeros((self.nmodes,15))
            bmean = empty((self.nmodes,2))
            hpol_avg = empty(self.nmodes)
            tn = self.dhpol
            ##----------------------------------------------------
            # Zhongwei: zzz component only
            if component == 'zzz':
               print('Considering only ZZZ component')
               for i in range(0,len(tn)):
                   for j in range(3):
                       for k in range(3):
                           for l in range(3):
                               if j == k == l == 2:
                                  pass
                               else:
                                  tn[i][j][k][l] = 0
            #----------------------------------------------------
            cj = tn.conjugate()
        else:
            b_avg = zeros((1,15))
            bmean = empty((1,2))
            hpol_avg = empty(1)
            table = ['STATIC', 'SHG', 'EOPE', 'OR', 'FD']
            for item in table:
                if item in self.hyperpolarizability.keys():
                    if real:
                      tn = array([self.hyperpolarizability[item].real],dtype=complex)
                    elif imag:
                      tn = array([self.hyperpolarizability[item].imag],dtype=complex)
                    else:
                      tn = array([self.hyperpolarizability[item]],dtype=complex)
                    cj = tn.conjugate()
                    break

        # Defining the coordinates
        first = range(3)
        second = range(3)
        third = range(3)

        for n in range(len(tn)):
            #print(n, self.v_frequencies[n])
            # b1
            b_avg[n][0] = sum([tn[n][i][i][i]*cj[n][i][i][i]
                               for i in first])
            # Components of the dominant bterm contribution (coeffcients involved) to both
            # <\beta_ZZZ> and <\beta_XZZ> for the mode at 1591.896 wavenumber^{-1} of CV
            if ratio and cv_mode and mode_bterm:
               if n == 96:
                  print ('')
                  print ('ZZZ')
                  print ('i, j, beta_aaa*cj(beta_aaa)')
                  for i in first:
                      print (i, (tn[n][i][i][i]*cj[n][i][i][i]).real*1/7.0)

               if n == 96:
                  print ('')
                  print ('XZZ')
                  print ('i, j, beta_aaa*cj(beta_aaa)')
                  for i in first:
                      print (i, (tn[n][i][i][i]*cj[n][i][i][i]).real*1/35.0)
            #b1 = 0.00
            #for i in first:
            #    b1 = b1 + tn[n][i][i][i]*cj[n][i][i][i]
            #print('b1', b_avg[n][0], b1)

            # b2
            b_avg[n][1] = sum([tn[n][i][i][j]*cj[n][i][i][j]
                               if i != j else 0.0 for i in first for j in second])
            #b2 = 0.00
            #for i in first:
            #    for j in second:
            #        if i == j:
            #            continue
            #        else:
            #            b2 = b2 + tn[n][i][i][j]*cj[n][i][i][j]
            #print('b2', b_avg[n][1], b2)

            # b3
            b_avg[n][2] = sum([tn[n][i][i][i]*cj[n][i][j][j]
                               if i != j else 0.0 for i in first for j in second])
            #b3 = 0.00
            #for i in first:
            #    for j in second:
            #        if i == j:
            #            continue
            #        else:
            #            b3 = b3 + tn[n][i][i][i]*cj[n][i][j][j]
            #print('b3', b_avg[n][2], b3)

            # b4
            b_avg[n][3] = sum([tn[n][j][i][i]*cj[n][i][i][j]
                               if i != j else 0.0 for i in first for j in second])
            #b4 = 0.00
            #for i in first:
            #    for j in second:
            #        if i == j:
            #            continue
            #        else:
            #            b4 = b4 + tn[n][j][i][i]*cj[n][i][i][j]
            #print('b4', b_avg[n][3], b4)

            # b5
            b_avg[n][4] = sum([tn[n][i][i][i]*cj[n][j][j][i]
                               if i != j else 0.0 for i in first for j in second])
            #b5 = 0.00
            #for i in first:
            #    for j in second:
            #        if i == j:
            #            continue
            #        else:
            #            b5 = b5 + tn[n][i][i][i]*cj[n][j][j][i]
            #print('b5', b_avg[n][4], b5)

            # b6
            b_avg[n][5] = sum([tn[n][j][i][i]*cj[n][j][i][i]
                               if i != j else 0.0 for i in first for j in second])
            #b6 = 0.00
            #for i in first:
            #    for j in second:
            #        if i == j:
            #            continue
            #        else:
            #            b6 = b6 + tn[n][j][i][i]*cj[n][j][i][i]
            #print('b6', b_avg[n][5], b6)

            # b7
            b_avg[n][6] = sum([tn[n][i][i][j]*cj[n][j][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #b7 = 0.00
            #for i in first:
            #    for j in second:
            #        for k in third:
            #            if (i == j or i == k or j == k):
            #                continue
            #            else:
            #                b7 = b7 + tn[n][i][i][j]*cj[n][j][k][k]
            #print('b7', b_avg[n][6], b7)

            # b8
            b_avg[n][7] = sum([tn[n][j][i][i]*cj[n][j][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #b8 = 0.00
            #for i in first:
            #    for j in second:
            #        for k in third:
            #            if (i == j or i == k or j == k):
            #                continue
            #            else:
            #                b8 = b8 + tn[n][j][i][i]*cj[n][j][k][k]
            #print('b8', b_avg[n][7], b8)

            # b9
            b_avg[n][8] = sum([tn[n][i][i][j]*cj[n][k][k][j]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #b9 = 0.00
            #for i in first:
            #    for j in second:
            #        for k in third:
            #            if (i == j or i == k or j == k):
            #                continue
            #            else:
            #                b9 = b9 + tn[n][i][i][j]*cj[n][k][k][j]
            #print('b9', b_avg[n][8], b9)

            # b10
            b_avg[n][9] = sum([tn[n][i][j][k]*cj[n][i][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #b10 = 0.00
            #for i in first:
            #    for j in second:
            #        for k in third:
            #            if (i == j or i == k or j == k):
            #                continue
            #            else:
            #                b10 = b10 + tn[n][i][j][k]*cj[n][i][j][k]
            #print('b10', b_avg[n][9], b10)

            # b11
            b_avg[n][10] = sum([tn[n][i][j][k]*cj[n][j][i][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #b11 = 0.00
            #for i in first:
            #    for j in second:
            #        for k in third:
            #            if (i == j or i == k or j == k):
            #                continue
            #            else:
            #                b11 = b11 + tn[n][i][j][k]*cj[n][j][i][k]
            #print('b11', b_avg[n][10], b11)

            # b12
            b_avg[n][11] = sum([tn[n][i][j][j]*cj[n][i][j][j]
                               if i != j else 0.0
                               for i in first for j in second])
            # Components of the dominant bterm contribution (coeffcients involved)
            # to <\beta_XZZ> for the mode at 1591.896 wavenumber^{-1} of CV
            if ratio and cv_mode and mode_bterm:
               if n == 96:
                  print ('')
                  print ('XZZ')
                  print ('i, j, beta_abb*cj(beta_abb)')
                  for i in first:
                      for j in second:
                          if i!=j:
                             print (i,j,(tn[n][i][j][j]*cj[n][i][j][j]).real*3/35.0)
            #b12 = 0.00
            #for i in first:
            #    for j in second:
            #        if i == j:
            #            continue
            #        else:
            #            b12 = b12 + tn[n][i][j][j]*cj[n][i][j][j]
            #print('b12', b_avg[n][11], b12)

            # b13
            b_avg[n][12] = sum([tn[n][i][i][j]*cj[n][j][i][i]
                               if i != j else 0.0
                               for i in first for j in second])
            #b13 = 0.00
            #for i in first:
            #    for j in second:
            #        if i == j:
            #            continue
            #        else:
            #            b13 = b13 + tn[n][i][i][j]*cj[n][j][i][i]
            #print('b13', b_avg[n][12], b13)

            # b14
            b_avg[n][13] = sum([tn[n][i][j][j]*cj[n][i][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #b14 = 0.00
            #for i in first:
            #    for j in second:
            #        for k in third:
            #            if (i == j or i == k or j == k):
            #                continue
            #            else:
            #                b14 = b14 + tn[n][i][j][j]*cj[n][i][k][k]
            #print('b14', b_avg[n][13], b14)

            # b15
            b_avg[n][14] = sum([tn[n][i][i][k]*cj[n][j][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #b15 = 0.00
            #for i in first:
            #    for j in second:
            #        for k in third:
            #            if (i == j or i == k or j == k):
            #                continue
            #            else:
            #                b15 = b15 + tn[n][i][i][k]*cj[n][j][j][k]
            #print('b15', b_avg[n][14], b15)

            # Calculate averages
            #bmean1 = 0.00
            #bmean1 = b1/7
            #bmean1 = bmean1 + 4*b2/35 + 2*b3/35 + 4*b4/35 + 4*b5/35 + b6/35
            #bmean1 = bmean1 + 4*b7/105 + b8/105 + 4*b9/105 + 2*b10/105 + 4*b11/105

            #bmean2 = 0.00
            #bmean2 = b1/35
            #bmean2 = bmean2 + 4*b3/105 - 4*b5/70 + 8*b2/105 + 3*b12/35 - 4*b13/70
            #bmean2 = bmean2 + b14/35 -4*b15/210 - 4*b7/210 + 2*b10/35 -4*b11/210

            bmean[n][0] = ( b_avg[n][0]/7 + 4*b_avg[n][1]/35 + 2*b_avg[n][2]/35
                          + 4*b_avg[n][3]/35 + 4*b_avg[n][4]/35 + b_avg[n][5]/35
                          + 4*b_avg[n][6]/105 + b_avg[n][7]/105 + 4*b_avg[n][8]/105
                          + 2*b_avg[n][9]/105 + 4*b_avg[n][10]/105)
            bmean[n][1] = ( b_avg[n][0]/35 + 4*b_avg[n][2]/105 - 4*b_avg[n][4]/70
                          + 8*b_avg[n][1]/105 + 3*b_avg[n][11]/35 - 4*b_avg[n][12]/70
                          + b_avg[n][13]/35 - 4*b_avg[n][14]/210 - 4*b_avg[n][6]/210
                          + 2*b_avg[n][9]/35 - 4*b_avg[n][10]/210)
            #print('bmean1', bmean[n][0], bmean1)
            #print('bmean2', bmean[n][1], bmean2)
            if ratio:
               hpol_avg[n] = bmean[n][1]/bmean[n][0]
            else:
               hpol_avg[n] = bmean[n][0] + bmean[n][1]
            # Depolarization ratios for hyper-Raman scattering (useful for
            # comparing with the Quinet paper)
            #print(self.v_frequencies[n],bmean[n][1]/bmean[n][0],hpol_avg[n])

            # Depol ratio with all bterm contributions (no coefficients
            # involved) for the mode at 1591.896 wavenumber^{-1} of CV
            if ratio and cv_mode:
                  if n == 96:
                     print ('')
                     print ('Mode HRS_Depol_Ratio')
                     print (self.v_frequencies[n], hpol_avg[n])
                     print ('')
                     print ('No. Values')
                     for i in range(0, 15): 
                         print (i,':', b_avg[n][i])

        return hpol_avg

    def shpol_average(self):
        '''Calculates the average of the second hyperpolarizability for second hyper-Raman
        scattering collected perpendicular to the incident light direction.'''

        from numpy import zeros, empty, array

        assert 'SECOND HYPERPOLARIZABILITY' in self.calctype, ('shpol_average(): '
                                 'No second hyperpolarizabilities were collected.')

        # Arrays storing the total second hyperpolarizability averages.
        if 'FREQUENCIES' in self.calctype:
            g_avg = zeros((self.nmodes,29))
            gmean = empty((self.nmodes,2))
            shpol_avg = empty(self.nmodes)
            tn = self.dshpol
            cj = tn.conjugate()
        else:
            g_avg = zeros((1,29))
            gmean = empty((1,2))
            shpol_avg = empty(1)
            table = ['STATIC', 'THG', 'EFISHG', 'EFIOR', 'OKE', 'IDRI', 'FD']
            for item in table:
                if item in self.shyperpolarizability.keys():
                    tn = array([self.shyperpolarizability[item]],dtype=complex)
                    cj = tn.conjugate()
                    break

        # Defining the coordinates
        first = range(3)
        second = range(3)
        third = range(3)

        for n in range(len(tn)):
            #print(n, self.v_frequencies[n])
            #---------------------------------------------------------------------
            # g1, aaaa*aaaa
            g_avg[n][0] = sum([tn[n][i][i][i][i]*cj[n][i][i][i][i]
                               for i in first])
            #---------------------------------------------------------------------

            #---------------------------------------------------------------------
            # g2, aaab*aaab
            g_avg[n][1] = sum([tn[n][i][i][i][j]*cj[n][i][i][i][j]
                               if i != j else 0.0 for i in first for j in second])

            # g3, aaaa*aabb
            g_avg[n][2] = sum([tn[n][i][i][i][i]*cj[n][i][i][j][j]
                               if i != j else 0.0 for i in first for j in second])

            # g4, aaab*baaa
            g_avg[n][3] = sum([tn[n][i][i][i][j]*cj[n][j][i][i][i]
                               if i != j else 0.0 for i in first for j in second])

            # g5, baaa*baaa
            g_avg[n][4] = sum([tn[n][j][i][i][i]*cj[n][j][i][i][i]
                               if i != j else 0.0 for i in first for j in second])

            # g6, aaaa*baab
            g_avg[n][5] = sum([tn[n][i][i][i][i]*cj[n][j][i][i][j]
                               if i != j else 0.0 for i in first for j in second])
            #---------------------------------------------------------------------

            #---------------------------------------------------------------------
            # g7, aabb*aabb
            g_avg[n][6] = sum([tn[n][i][i][j][j]*cj[n][i][i][j][j]
                               if i != j else 0.0 for i in first for j in second])

            # g8, aaab*abbb
            g_avg[n][7] = sum([tn[n][i][i][i][j]*cj[n][i][j][j][j]
                               if i != j else 0.0 for i in first for j in second])

            # g9, abbb*baaa
            g_avg[n][8] = sum([tn[n][i][j][j][j]*cj[n][j][i][i][i]
                               if i != j else 0.0 for i in first for j in second])

            # g10, aabb*baab
            g_avg[n][9] = sum([tn[n][i][i][j][j]*cj[n][j][i][i][j]
                               if i != j else 0.0 for i in first for j in second])

            # g11, baab*baab
            g_avg[n][10] = sum([tn[n][j][i][i][j]*cj[n][j][i][i][j]
                               if i != j else 0.0 for i in first for j in second])

            # g12, baaa*babb (one possible typo says: "aaab*bbba")
            g_avg[n][11] = sum([tn[n][j][i][i][i]*cj[n][j][i][j][j]
                               if i != j else 0.0 for i in first for j in second])

            # g13, aaab*bbba (one possible typo says: "baaa*bbba")
            g_avg[n][12] = sum([tn[n][i][i][i][j]*cj[n][j][j][j][i]
                               if i != j else 0.0 for i in first for j in second])

            # g14, aaaa*bbbb
            g_avg[n][13] = sum([tn[n][i][i][i][i]*cj[n][j][j][j][j]
                               if i != j else 0.0 for i in first for j in second])
            #---------------------------------------------------------------------

            #---------------------------------------------------------------------
            # g15, aabc*aabc
            g_avg[n][14] = sum([tn[n][i][i][j][k]*cj[n][i][i][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g16, aabb*aacc
            g_avg[n][15] = sum([tn[n][i][i][j][j]*cj[n][i][i][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g17, aaac*abbc
            g_avg[n][16] = sum([tn[n][i][i][i][k]*cj[n][i][j][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g18, abcc*baaa
            g_avg[n][17] = sum([tn[n][i][j][k][k]*cj[n][j][i][i][i]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g19, aacc*baab
            g_avg[n][18] = sum([tn[n][i][i][k][k]*cj[n][j][i][i][j]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g20, aabc*baac
            g_avg[n][19] = sum([tn[n][i][i][j][k]*cj[n][j][i][i][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g21, baac*baac
            g_avg[n][20] = sum([tn[n][j][i][i][k]*cj[n][j][i][i][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g22, aaac*babc
            g_avg[n][21] = sum([tn[n][i][i][i][k]*cj[n][j][i][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g23, aaab*bacc
            g_avg[n][22] = sum([tn[n][i][i][i][j]*cj[n][j][i][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g24, baaa*bacc
            g_avg[n][23] = sum([tn[n][j][i][i][i]*cj[n][j][i][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g25, baaa*bacc
            g_avg[n][24] = sum([tn[n][i][i][i][i]*cj[n][j][j][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g26, babc*caaa
            g_avg[n][25] = sum([tn[n][j][i][j][k]*cj[n][k][i][i][i]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g27, baac*caab
            g_avg[n][26] = sum([tn[n][j][i][i][k]*cj[n][k][i][i][j]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g28, baab*caac
            g_avg[n][27] = sum([tn[n][j][i][i][j]*cj[n][k][i][i][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

            # g29, baaa*cabc (one possible typo says: "aaaa*cbbc")
            g_avg[n][28] = sum([tn[n][j][i][i][i]*cj[n][k][i][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])
            #---------------------------------------------------------------------

            # Calculate averages
            #gmean1 = g1/9
            #gmean1 = gmean1 + g2/7 + 2*g3/21 + 2*g4/21 + g5/63 + 2*g6/21
            #gmean1 = gmean1 + 3*g7/35 + 2*g8/35 + 2*g9/105 + 6*g10/35 + 3*g11/35 
            #                + 2*g12/35 + 6*g13/35 + 2*g14/105 
            #gmean1 = gmean1 + 2*g15/35 + g16/35 + 2*g17/35 + 2*g18/105 + 2*g19/35
            #                + 4*g20/35 + g21/105 + 4*g22/35 + 2*g23/35 + 2*g24/105
            #                + g25/63 + 4*g26/105 + g27/35 + g28/35 + g29/315

            #gmean2 = g1/63
            #gmean2 = gmean2 + 2*g2/35 + 4*g3/105 - g4/21 + 4*g5/63 - g6/21
            #gmean2 = gmean2 + 3*g7/35 + 2*g8/35 - g9/105 - 3*g10/35 + 3*g11/35
            #                - 3*g12/35 + 2*g13/35 - g14/105 
            #gmean2 = gmean2 + 2*g15/35 + g16/35 + 2*g17/35 - g18/105 - g19/35
            #                - 2*g20/35 + 4*g21/105 (a possible type says"4*g21/35)
            #                - 2*g22/35 -g23/35 + 8*g24/105 - g25/105 - g26/210
            #                - g27/70 - g28/70 - g29/70

            gmean[n][0] = ( 1*g_avg[n][0]/9    + 1*g_avg[n][1]/7    + 2*g_avg[n][2]/21
                          + 2*g_avg[n][3]/21   + 1*g_avg[n][4]/63   + 2*g_avg[n][5]/21
                          + 3*g_avg[n][6]/35   + 2*g_avg[n][7]/35   + 2*g_avg[n][8]/105
                          + 6*g_avg[n][9]/35   + 3*g_avg[n][10]/35  + 2*g_avg[n][11]/35
                          + 6*g_avg[n][12]/35  + 2*g_avg[n][13]/105 + 2*g_avg[n][14]/35
                          + 1*g_avg[n][15]/35  + 2*g_avg[n][16]/35  + 2*g_avg[n][17]/105
                          + 2*g_avg[n][18]/35  + 4*g_avg[n][19]/35  + 1*g_avg[n][20]/105
                          + 4*g_avg[n][21]/35  + 2*g_avg[n][22]/35  + 2*g_avg[n][23]/105
                          + 1*g_avg[n][24]/63  + 4*g_avg[n][25]/105 + 1*g_avg[n][26]/35
                          + 1*g_avg[n][27]/35  + 1*g_avg[n][28]/315 )

            gmean[n][1] = ( 1*g_avg[n][0]/63   + 2*g_avg[n][1]/35   + 4*g_avg[n][2]/105
                          - 1*g_avg[n][3]/21   + 4*g_avg[n][4]/63   - 1*g_avg[n][5]/21
                          + 3*g_avg[n][6]/35   + 2*g_avg[n][7]/35   - 1*g_avg[n][8]/105
                          - 3*g_avg[n][9]/35   + 3*g_avg[n][10]/35  - 3*g_avg[n][11]/35
                          + 2*g_avg[n][12]/35  - 1*g_avg[n][13]/105 + 2*g_avg[n][14]/35
                          + 1*g_avg[n][15]/35  + 2*g_avg[n][16]/35  - 1*g_avg[n][17]/105
                          - 1*g_avg[n][18]/35  - 2*g_avg[n][19]/35  + 4*g_avg[n][20]/105
                          - 2*g_avg[n][21]/35  - 1*g_avg[n][22]/35  + 8*g_avg[n][23]/105
                          - 1*g_avg[n][24]/105 + 1*g_avg[n][25]/210 - 1*g_avg[n][26]/70
                          - 1*g_avg[n][27]/70  - 1*g_avg[n][28]/70 )

            #print('gmean1', gmean[n][0], gmean1)
            #print('gmean2', gmean[n][1], gmean2)
            shpol_avg[n] = gmean[n][0] + gmean[n][1]
            # Depolarization ratios for second hyper-Raman scattering 
            #print(self.v_frequencies[n],gmean[n][1]/gmean[n][0],shpol_avg[n])

        return shpol_avg

    def pol_minmax(self, property=None, dim=False, qm=True):
        '''Returns the minimum and maximum isotropic polarizabilities
        of a group of polarizabilities.

        The option *property* can be either 'pol' or 'ord' for polarizability
        or optical rotation, respectively.  If neither is given, the property
        one will be guessed based on the calctype.

        If the boolean *dim* is true, then the min and max of the DIM
        system will be returned.  Same for the boolean *qm*.
        If both are true, the polarizabilities of the two systems are summed.

        Raises :py:exc:`AssertionError` if no polarizabilities were collected.

        '''
        return (self.isotropic(property, qm=qm, dim=dim).min(),
                self.isotropic(property, qm=qm, dim=dim).max())


    def pol_diagonalize(self, property=None, dim=False):
        '''Diagonalizes the :py:attr:`polarizability tensor <polarizability>`
        in place.  The rotation matrix after diagonalization is returned.

        The option *property* can be either 'pol' or 'ord' for polarizability
        or optical rotation, respectively.  If neither is given, the property
        one will be guessed based on the calctype.

        If the boolean *dim* is True, then
        :py:attr:`DIM polarizability <dim_pol>` will be diagonalized INSTEAD
        of the qm polarizability.

        Raises :py:exc:`AssertionError` if no polarizabilities were collected.

        '''
        from numpy import diagflat, zeros_like
        from numpy.linalg import eig, LinAlgError

        if property == 'pol':
            assert 'POLARIZABILITY' in self.calctype, ('pol_diagonalize(): '
                                         'No polarizabilities were collected.')
        elif property == 'ord':
            assert 'OPTICAL ROTATION' in self.calctype, ('pol_diagonalize(): '
                                        'No optical rotations were collected.')
        elif property is None:
            # Try to guess the property
            if 'POLARIZABILITY' in self.calctype:
                property = 'pol'
            elif 'OPTICAL ROTATION' in self.calctype:
                property = 'ord'
            else:
                raise AssertionError ('pol_diagonalize(): Nothing to do.')
        if dim: assert 'DIM' in self.calctype, ('pol_diagonalize(): '
                                                    'Not a DIMQM calculation.')

        if dim:
            tn = self.dim_pol
        else:
            if property == 'pol':
                tn = self.polarizability
            elif property == 'ord':
                tn = self.ord

        # Find the eigenvector and eigenvalues of the tensor
        # The eigenvector is the diagonal of the tensor, and the
        # eigenvalues are the rotation matrix.
        if rot:
           rotation = zeros_like(tn)

        for n in range(self.npol):
            try:
                eigenvec, eigenval = eig(tn[n])
            except LinAlgError:
                raise ChemDataError (
                  'Diagonalization of polarizability tensor does not converge')
            else:
                tn[n] = diagflat(eigenvec)
                if rot:
                    rotation[n] = eigenval

        # Remember the state
        self._diagonalized = True

        return rotation

    def fold_Atensor(self, action=None, astensor=False):
        '''Folds the A-tensor from a 3x6 to a 3x3x3 tensor, or vice versa.
        The keyword :py:attr:`action` can be either :py:attr:`FOLD`, which
        folds the 3x3 quadrupole into a length 6 quadrupole representation,
        or :py:attr:`UNFOLD`, which unfolds the length 6 quadrupole operator
        into a 3x3 operator.'''
        from numpy import zeros
        if ((astensor and self.astensor is not None) or
            (not astensor and  self.atensor is not None)):
            if astensor:
                dtype = self.astensor.dtype
                ndim = len(self.astensor[0][0])
                A = self.astensor[:]
            else:
                dtype = self.atensor.dtype
                ndim = len(self.atensor[0][0])
                A = self.atensor[:]
            if action is None:
                if ndim == 6: action = 'UNFOLD'
                if ndim == 3: action = 'FOLD'
            action = action.upper()
            if action == 'UNFOLD' and ndim == 6:
                temp = zeros((len(A),3,3,3), dtype=dtype)
            elif action == 'FOLD' and ndim == 3:
                temp = zeros((len(A),3,6), dtype=dtype)
            else:
                return
            for mode in range(len(A)):
              for a in range(3):
                for i in range(6):
                  j=0; k=i
                  if i==3: j=1; k=1
                  if i==4: j=1; k=2
                  if i==5: j=2; k=2
                  if action == 'UNFOLD' and ndim == 6:
                    temp[mode][a][j][k] = A[mode][a][i]
                    temp[mode][a][k][j] = A[mode][a][i]
                  if action == 'FOLD' and ndim == 3:
                    temp[mode][a][i] = ( A[mode][a][j][k]
                                       + A[mode][a][k][j] ) / 2.
            try:
                if astensor:
                    self.astensor = temp
                else:
                    self.atensor = temp
            except UnboundLocalError:
                pass

    def fold_Dtensor(self, action=None, dstensor=False):
        '''Folds the D-tensor from a 6x3 to a 3x3x3 tensor, or vice versa.
        The keyword :py:attr:`action` can be either :py:attr:`FOLD`, which
        folds the 3x3 quadrupole operator into a length 6 quadrupole representation,
        or :py:attr:`UNFOLD`, which unfolds the length 6 quadrupole operator
        into a 3x3 operator.'''
        from numpy import zeros
        if ((dstensor and self.dstensor is not None) or
            (not dstensor and self.dtensor is not None)):
            if dstensor:
                dtype = self.dstensor.dtype
                ndim = len(self.dstensor[0][0])
                D = self.dstensor[:]
            else:
                dtype = self.dtensor.dtype
                ndim = len(self.dtensor[0])
                D = self.dtensor[:]
            if action is None:
                if ndim == 6: action = 'UNFOLD'
                if ndim == 3: action = 'FOLD'
            action = action.upper()
            if action == 'UNFOLD' and ndim == 6:
                temp = zeros((len(D),3,3,3), dtype=dtype)
            elif action == 'FOLD' and ndim == 3:
                if dstensor:
                    temp = zeros((len(D),3,6), dtype=dtype)
                else:
                    temp = zeros((len(D),6,3), dtype=dtype)
            else:
                return
            for mode in range(len(D)):
              for a in range(3):
                for i in range(6):
                  j=0; k=i
                  if i==3: j=1; k=1
                  if i==4: j=1; k=2
                  if i==5: j=2; k=2
                  if action == 'UNFOLD' and ndim == 6:
                    if dstensor:
                        temp[mode][a][j][k] = D[mode][a][i]
                        temp[mode][a][k][j] = D[mode][a][i]
                    else:
                        temp[mode][j][k][a] = D[mode][i][a]
                        temp[mode][k][j][a] = D[mode][i][a]
                  if action == 'FOLD' and ndim == 3:
                    if dstensor:
                        temp[mode][a][i] = ( D[mode][a][j][k]
                                           + D[mode][a][k][j] ) / 2.
                    else:
                        temp[mode][i][a] = ( D[mode][j][k][a]
                                           + D[mode][k][j][a] ) / 2.
            try:
                if dstensor:
                    self.dstensor = temp
                else:
                    self.dtensor = temp
            except UnboundLocalError:
                pass

    def fold_Ctensor(self, action=None):
        '''Folds the Ctensor from a 6x6 to a 3x3x3x3 tensor, or vice versa.
        The keyword :py:attr:`action` can be either :py:attr:`FOLD`, which
        folds the 3x3 quadrupole operators into a length 6 quadrupole representation,
        or :py:attr:`UNFOLD`, which unfolds the length 6 quadrupole operators
        into 3x3 operators.'''

        from numpy import zeros
        if self.ctensor is not None:
            dtype = self.ctensor.dtype
            ndim = len(self.ctensor[0])
            if action is None:
                if ndim == 6: action = 'UNFOLD'
                if ndim == 3: action = 'FOLD'
            action = action.upper()
            if action == 'UNFOLD' and ndim == 6:
                temp = zeros((len(self.ctensor),3,3,3,3), dtype=dtype)
            elif action == 'FOLD' and ndim == 3:
                temp = zeros((len(self.ctensor),6,6), dtype=dtype)
            else:
                return
            for mode in range(len(self.ctensor)):
              for a in range(6):
                b=0; c=a
                if a==3: b=1; c=1
                if a==4: b=1; c=2
                if a==5: b=2; c=2
                for i in range(6):
                  j=0; k=i
                  if i==3: j=1; k=1
                  if i==4: j=1; k=2
                  if i==5: j=2; k=2
                  if action == 'UNFOLD' and ndim == 6:
                    temp[mode][b][c][j][k] = self.ctensor[mode][a][i]
                    temp[mode][b][c][k][j] = self.ctensor[mode][a][i]
                    temp[mode][c][b][j][k] = self.ctensor[mode][a][i]
                    temp[mode][c][b][k][j] = self.ctensor[mode][a][i]
                  elif action == 'FOLD' and ndim == 3:
                    temp[mode][a][i] = ( self.ctensor[mode][b][c][j][k]
                                       + self.ctensor[mode][b][c][k][j]
                                       + self.ctensor[mode][c][b][j][k]
                                       + self.ctensor[mode][c][b][k][j] ) / 4.
            try:
                self.ctensor = temp
            except UnboundLocalError:
                pass


    def fold_Btensor(self, action=None):
        '''Folds (unfolds) the B-tensor from a 3x10 to a 3x3x3x3 tensor (or vice versa).
        The keyword :py:attr:`action` can be either :py:attr:`FOLD`, which
        folds the 3x3x3 octupole operator into a length 10 octupole representation,
        or :py:attr:`UNFOLD`, which unfolds the length 10 octupole operator
        into a 3x3x3 operator.'''

        from numpy import zeros
        if self.btensor is not None:
            dtype = self.btensor.dtype
            ndim = len(self.btensor[0][0])
            if action is None:
                if ndim == 10: action = 'UNFOLD'
                if ndim == 3: action = 'FOLD'
            action = action.upper()
            if action == 'UNFOLD' and ndim == 10:
                temp = zeros((len(self.btensor),3,3,3,3), dtype=dtype)
            elif action == 'FOLD' and ndim == 3:
                temp = zeros((len(self.btensor),3,10), dtype=dtype)
            else:
                return
            for mode in range(len(self.btensor)):
              for i in range(3):
                for j in range(10):
                  if j==0: a=0; b=0; c=0
                  if j==1: a=0; b=0; c=1
                  if j==2: a=0; b=0; c=2
                  if j==3: a=0; b=1; c=1
                  if j==4: a=0; b=1; c=2
                  if j==5: a=0; b=2; c=2
                  if j==6: a=1; b=1; c=1
                  if j==7: a=1; b=1; c=2
                  if j==8: a=1; b=2; c=2
                  if j==9: a=2; b=2; c=2
                  if action == 'UNFOLD' and ndim == 10:
                    temp[mode][i][a][b][c] = self.btensor[mode][i][j]
                    temp[mode][i][a][c][b] = self.btensor[mode][i][j]
                    temp[mode][i][b][a][c] = self.btensor[mode][i][j]
                    temp[mode][i][b][c][a] = self.btensor[mode][i][j]
                    temp[mode][i][c][a][b] = self.btensor[mode][i][j]
                    temp[mode][i][c][b][a] = self.btensor[mode][i][j]
                  elif action == 'FOLD' and ndim == 3:
                    temp[mode][i][j] = ( self.btensor[mode][i][a][b][c]
                                       + self.btensor[mode][i][a][c][b]
                                       + self.btensor[mode][i][b][a][c]
                                       + self.btensor[mode][i][b][c][a]
                                       + self.btensor[mode][i][c][a][b]
                                       + self.btensor[mode][i][c][b][a] ) / 6.
            try:
                self.btensor = temp
            except UnboundLocalError:
                pass


    def rotate_polarizabilities(self, R):
        '''Rotate the higher order polarizabilities such as
        the A, B, C, D, and G tensors.
        *R* is the rotation matrix.'''

        from numpy import einsum, array

        # Indices arrays for folding/unfolding tensors
        iQ = [[0, 1, 2], [1, 3, 4], [2, 4, 5]]
        iQQ = [tuple([0,0]), tuple([0,1]), tuple([0,2]), tuple([1,1]), tuple([1,2]), tuple([2,2])]

        # Rotate each tensor if present
        if self.qm_pol is not None:
            self.qm_pol = einsum('ij,kl,ajl', R, R, self.qm_pol)
        if self.ord is not None:
            self.ord = einsum('ij,kl,ajl', R, R, self.ord)
        if self.gtensor is not None:
            self.gtensor = einsum('ij,kl,ajl', R, R, self.gtensor)
        if self.atensor is not None:
            self.fold_Atensor('UNFOLD')
            self.atensor = einsum('ij,kl,mn,ajln', R, R, R, self.atensor)
        if self.dtensor is not None:
            self.fold_Dtensor('UNFOLD')
            self.dtensor = einsum('ij,kl,mn,ajln', R, R, R, self.dtensor)
        if self.ctensor is not None:
            self.fold_Ctensor('UNFOLD')
            self.ctensor = einsum('ij,kl,mn,op,ajlnp', R, R, R, R, self.ctensor)
        if self.btensor is not None:
            self.fold_Btensor('UNFOLD')
            self.btensor = einsum('ij,kl,mn,op,ajlnp', R, R, R, R, self.btensor)
        if self.hyperpolarizability is not None:
            self.hyperpolarizability = einsum('ia,jb,kc,zabc->zijk', R, R, R, self.hyperpolarizability)
#### implementing rotation of hyperpol now
        if self.dhpol is not None:
            self.dhpol = einsum('ia,jb,kc,zabc->zijk', R, R, R, self.dhpol)
##############################################
        if self.beta_ddq is not None:
            temp = array([[[[[self.beta_ddq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in range(3)]
                         for j in range(3)] for i in range(self.nmodes)])
            temp = einsum('ia,jb,kc,ld,zabcd->zijkl', R, R, R, R, temp)
            self.beta_ddq = array([[[[temp[i][j][k][l] for l in iQQ] for k in range(3)] for j in range(3)]
                                  for i in range(self.nmodes)])
        if self.beta_dqd is not None:
            temp = array([[[[[self.beta_dqd[i][j][k][l] for l in range(3)] for k in iQ[m]] for m in range(3)]
                         for j in range(3)] for i in range(self.nmodes)])
            temp = einsum('ia,jb,kc,ld,zabcd->zijkl', R, R, R, R, temp)
            self.beta_dqd = array([[[[temp[i][j][k][l] for l in range(3)] for k in iQQ] for j in range(3)]
                                  for i in range(self.nmodes)])
        if self.beta_qdd is not None:
            temp = array([[[[[self.beta_qdd[i][j][k][l] for l in range(3)] for k in range(3)] for j in iQ[m]]
                         for m in range(3)] for i in range(self.nmodes)])
            temp = einsum('ia,jb,kc,ld,zabcd->zijkl', R, R, R, R, temp)
            self.beta_qdd = array([[[[temp[i][j][k][l] for l in range(3)] for k in range(3)] for j in iQQ]
                                  for i in range(self.nmodes)])
        if self.beta_dqq is not None:
            temp = array([[[[[[self.beta_dqq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in iQ[n]]
                         for n in range(3)] for j in range(3)] for i in range(self.nmodes)])
            temp = einsum('ia,jb,kc,ld,me,zabcde->zijklm', R, R, R, R, R, temp)
            self.beta_dqq = array([[[[temp[i][j][k][l] for l in iQQ] for k in iQQ] for j in range(3)]
                                  for i in range(self.nmodes)])
        if self.beta_qdq is not None:
            temp = array([[[[[[self.beta_qdq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in range(3)]
                         for j in iQ[n]] for n in range(3)] for i in range(self.nmodes)])
            temp = einsum('ia,jb,kc,ld,me,zabcde->zijklm', R, R, R, R, R, temp)
            self.beta_qdq = array([[[[temp[i][j][k][l] for l in iQQ] for k in range(3)] for j in iQQ]
                                  for i in range(self.nmodes)])
        if self.beta_qqd is not None:
            temp = array([[[[[[self.beta_qqd[i][j][k][l] for l in range(3)] for k in iQ[m]] for m in range(3)]
                         for j in iQ[n]] for n in range(3)] for i in range(self.nmodes)])
            temp = einsum('ia,jb,kc,ld,me,zabcde->zijklm', R, R, R, R, R, temp)
            self.beta_qqd = array([[[[temp[i][j][k][l] for l in range(3)] for k in iQQ] for j in iQQ]
                                  for i in range(self.nmodes)])
        if self.beta_qqq is not None:
            temp = array([[[[[[[self.beta_qqq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in iQ[n]]
                         for n in range(3)] for j in iQ[o]] for o in range(3)] for i in range(self.nmodes)])
            temp = einsum('ia,jb,kc,ld,me,nf,zabcdef->zijklmn', R, R, R, R, R, R, temp)
            self.beta_qqq = array([[[[temp[i][j][k][l] for l in iQQ] for k in iQQ] for j in iQQ]
                                  for i in range(self.nmodes)])
