from __future__ import print_function, division
from .errorclass import ChemDataError
import sys

class Raman_IR(object):
    '''Extends the ChemData class with methods to manipulate Raman and IR.'''

    def printVibSpec(self, laser=514.5):
        '''Pretty-print the vibrational spectroscopic data to
        standard output.

        '''
        from numpy import vstack
        from numpy import fastCopyAndTranspose as fcat

        assert 'FREQUENCIES' in self.calctype, ('printVibSpec(): '
                                              'Not a FREQUENCIES calculation.')

        # Check if this will be IR only, or Raman and IR
        if 'RAMAN' in self.calctype:
            head = ' {0}  {1}  {2}  {3}'.format('Frequency',
                                               'Cross Section',
                                               'Scattering Factor',
                                               'IR Coeff')
            f = ' {0[0]:9.3f}  {0[1]:<14.4E} {0[2]:<17.6g}  {0[3]:<7.5g}'
            data = fcat(vstack((self.v_frequencies,
                                self.cross_section(property='raman'),
                                self.scattering_factor(),
                                self.IR)))
        else:
            head = ' {0}  {1}'.format('Frequency', 'IR Coeff')
            f = ' {0[0]:9.3f}  {0[1]:<7.5g}'
            data = fcat(vstack((self.v_frequencies, self.IR)))

        # Now print
        print()
        print(head)
        for i in range(self.nmodes):
            print(f.format(data[i]))
        print()

    #def collect_raman_derivatives(self, dir=None, excitation=None, poop=False,
    #                              sR=0.01, gdipder=False, tdipder=False,
    #                              tdipref=None, hpol=False, shpol=False):
    # Zhongwei, provide the option to choose the 1st or 2nd or 3rd ... group of
    # polarizability tensors from higher-order property calculations, default
    # is the last group, i.e., "group = -1".
    def collect_raman_derivatives(self, dir=None, excitation=None, poop=False,
                                  sR=0.01, gdipder=False, tdipder=False,
                                  tdipref=None, hpol=False, shpol=False,
                                  group=-1, **kwargs):
        '''Collects all the :py:attr:`polarizability` or excitation derivative
        calculation files in a directory.

        The directory is assumed to be the current directory if none is given.
        The name of the files must be mode####.##-[mp].out.

        The option *excitation* chooses an excitation to take the derivative
        of.  The information will be put into :py:attr:`deltas`.  No
        polarizabilities will be collected.  The excitation is given as
        'n sym', which is the nth excitation of symmetry
        :py:attr:`sym <excitation_symmetries>`.

        If the option *excitation* is omitted, then the
        :py:attr:`polarizability` derivatives will be calculated.

        The self object is assumed to contain information collected from a
        frequency calculation (i.e. is of calctype **FREQUENCIES**), and then
        the derivative information is added to that base of information.
        This will add **RAMAN** to the calctype.

        Note that this method will not collect the DIM system values, even
        if they exist, because they will not change with normal modes.

        Set *poop* to True if you are using files that were generated before
        Dan fixed the stepsize issue.

        The option *sR*, which is the unweighted step size that was
        used to make the plus and minus direction coordinates for the
        3-point numerical differentiation, is defaulted to 0.01.  This
        value is used to calculate the mass-weighted step size.

        The option *gdipder* is used to calculate derivatives of the ground
        state dipole moment.

        The option *tdipder* is used to calculate derivatives of the transition
        dipole moment.

        The option *tdipref* defines the signs of the transition dipole moment
        from a reference calculation.  These are required because the sign
        of this quantity depends on the phase of the ground state wavefunction.

        The option *hpol* defines that hyperpolarizability derivatives are being
        calculated.

        The option *shpol* defines that second hyperpolarizability derivatives are
        being calculated.

        Raises :py:exc:`AssertionError` if this is not a
        **FREQUENCIES** calculation.
        '''
        from glob import glob
        from . import collect
        from natsort import natsort_key
        from numpy import append, array, where, ones, reshape, absolute, argmax
        from numpy import round as rnd
        import os
        # Verify that we aren't doing something dumb.
        assert 'FREQUENCIES' in self.calctype, ('collect_raman_derivatives(): '
                         'File', self.filename, "is not of type 'FREQUENCIES'")

        # Grab the files in the directory of interest.
        if dir is None:
            dir = os.curdir
        pfiles = glob(os.path.join(dir, 'mode*-p.out'))
        mfiles = glob(os.path.join(dir, 'mode*-m.out'))

        # Make sure that the number of plus and minus files is the same
        assert len(mfiles) == len(pfiles), ('collect_raman_derivatives(): '
                    'Number of minus and plus direction files does not match.')

        # Sort the filenames 'naturally', i.e. in numerical order, not ASCII.
        pfiles.sort(key=natsort_key)
        mfiles.sort(key=natsort_key)

        # Make a list of just the numbers (as strings)
        nums = []
        for file in pfiles:
            name = os.path.split(file)[1]
            num = name[4:] # remove 'mode' from front
            nums.append(num.split('-')[0]) # Keep only the number part

        # Loop over the numbers and collect the info into ChemData objects
        # Replace or augment some of the data.
        vfreq = array([], dtype=float)
        self.e_frequencies = array([], dtype=float)
        self.qm_pol = None
        self.dim_pol = None
        self.deltas = array([], dtype=float)
        self.dgdip = array([], dtype=float)
        self.dtdip = array([], dtype=float)
        self.TQM = array([], dtype=float)
        self.MDM = array([], dtype=float)
        self.dhpol = None
        self.dshpol = None
        self.b_e_frequencies = array([], dtype=float)
        for num in nums:
            pname = 'mode' + '-'.join((num, 'p')) + '.out'
            mname = 'mode' + '-'.join((num, 'm')) + '.out'
            pfile = os.path.join(dir, pname)
            mfile = os.path.join(dir, mname)
            try:
#------------------------- jbecca START --------------------------------
#               p = collect(pfile,project='SERS')
#               m = collect(mfile,project='SERS')
#       reverted because this routine is also used for hpol, 
#       which was broken thanks to this. Can be reimplemented
#       when better thought out
                p = collect(pfile) 
                m = collect(mfile)
                # For p/m files, SERS only needs a few attributes, not everything. --Pengchong, Nov. 2016                
#------------------------ jbecca END -----------------------------------
                m = collect(mfile)
            except IOError as e:
                sys.exit(str(e))

            if excitation:
                i, sym = excitation.split()
                try:
                    i = int(i) - 1
                except TypeError:
                    sys.exit('Requested excitation format invalid '+excitation)
                try:
                    ip = where(p.excitation_symmetries == sym)[0][i]
                    im = where(m.excitation_symmetries == sym)[0][i]
                except IndexError:
                    sys.exit('Requested excitation not found: '+excitation)
                temp = p.excitation_energies[ip] - m.excitation_energies[im]
                # Add this value to the array
                self.deltas = append(self.deltas, temp)
                # Ground state dipole moment derivatives
                if gdipder:
                    temp = p.dipole - m.dipole
                    self.dgdip = append(self.dgdip, temp)
                # Transition dipole moment derivatives
                # Also collect the transition quadrupole
                # and transition magnetic dipole moment
                # derivatives at the same time
                if tdipder:
                    # Make the sign of transition dipole moment from
                    # the plus file consistent with the equilibrium
                    # geometry transition dipole moment.
                    temp = absolute(p.TDM[ip])
                    temp = temp.argmax(axis=0)
                    #gtd = collect("tddft.out")
                    #sign = gtd.TDM[ip][temp]/absolute(gtd.TDM[ip][temp])
                    temp = p.TDM[ip][temp]/absolute(p.TDM[ip][temp])
                    sign = tdipref[ip]
                    if temp != sign:
                        p.TDM[ip] = -p.TDM[ip]
                        if p.TQM is not None: p.TQM[ip] = -p.TQM[ip]
                        if p.MDM is not None: p.MDM[ip] = -p.MDM[ip]
                    # Make the sign of transition dipole moment from
                    # the minus file consistent with the equilibrium
                    # geometry transition dipole moment.
                    temp = absolute(m.TDM[im])
                    temp = temp.argmax(axis=0)
                    #gtd = collect("tddft.out")
                    #sign = gtd.TDM[im][temp]/absolute(gtd.TDM[im][temp])
                    temp = m.TDM[im][temp]/absolute(m.TDM[im][temp])
                    sign = tdipref[im]
                    if temp != sign:
                        m.TDM[im] = -m.TDM[im]
                        if m.TQM is not None: m.TQM[im] = -m.TQM[im]
                        if m.MDM is not None: m.MDM[im] = -m.MDM[im]
                    if (p.TQM is not None and m.TQM is not None):
                        temp = p.TQM - m.TQM
                        self.TQM = append(self.TQM, temp)
                    if (p.MDM is not None and m.MDM is not None):
                        temp = p.MDM - m.MDM
                        self.MDM = append(self.MDM, temp)
                    temp = p.TDM - m.TDM
                    self.dtdip = append(self.dtdip, temp)
                    numexcite = p.nexcite
            elif hpol:
                # For static or frequency-dependent hyperpolarizabilities
                # (no damping).
                if ('STATIC' in p.calctype) or ('.SHG' in p.calctype):
                    try:
                        temp = array( p.hyperpolarizability['STATIC']
                                    - m.hyperpolarizability['STATIC'])
#jbecca: trying following block for dalton
#                    except KeyError:
#                        temp = array( p.hyperpolarizability['FD']
#                                    - m.hyperpolarizability['FD'])
                    except KeyError:
                        temp = array( p.hyperpolarizability['SHG']
                                    - m.hyperpolarizability['SHG'])
                    # Store hyperpolarizability derivatives, in units of
                    # a.u. / Angstrom**2 * amu**{1/2}
                    try:
                        self.dhpol = append(self.dhpol, [temp], axis=0)
                    except ValueError:
                        self.dhpol = [temp]
                # For damped frequency-dependent hyperpolarizabilities.
                # Although the derivatives are complex, this should not
                # be an issue since the code for generating hyper-Raman
                # spectra already makes use of complex conjugation.
                else:
                    ##----------------------------------------------------
                    ## Zhongwei: zzz component only
                    #for i in range(3):
                    #    for j in range(3):
                    #        for k in range(3):
                    #            if i == j == k == 2:
                    #               pass
                    #            else:
                    #               p.hyperpolarizability['FD'][i,j,k] = 0
                    #               m.hyperpolarizability['FD'][i,j,k] = 0
                    ##----------------------------------------------------
                    temp = array( p.hyperpolarizability['FD']
                                - m.hyperpolarizability['FD'])
                    # Store hyperpolarizability derivatives, in units of
                    # a.u. / Angstrom**2 * amu**{1/2}
                    try:
                        self.dhpol = append(self.dhpol, [temp], axis=0)
                    except ValueError:
                        self.dhpol = [temp]
                    # Add the laser wavelength
                    self.b_e_frequencies = append(self.b_e_frequencies,
                                                  p.b_e_frequencies[0])
            elif shpol:
                # For static or frequency-dependent second hyperpolarizabilities
                # (no damping).
                if 'STATIC' in p.calctype or 'THG' in p.calctype:
                    try:
                        temp = array( p.secondhyperpolarizability['STATIC']
                                    - m.secondhyperpolarizability['STATIC'])
                    except KeyError:
                        temp = array( p.secondhyperpolarizability['THG']
                                    - m.secondhyperpolarizability['THG'])
                    # Store second hyperpolarizability derivatives, in units of
                    # a.u. / Angstrom**2 * amu**{1/2} (to be edited)
                    try:
                        self.dshpol = append(self.dshpol, [temp], axis=0)
                    except ValueError:
                        self.dshpol = [temp]
                # For damped frequency-dependent second hyperpolarizabilities.
                # Although the derivatives are complex, this should not
                # be an issue since the code for generating second hyper-Raman
                # spectra already makes use of complex conjugation.
                else:
                    temp = array( p.secondhyperpolarizability['FD']
                                - m.secondhyperpolarizability['FD'])
                    # Store second hyperpolarizability derivatives, in units of
                    # a.u. / Angstrom**2 * amu**{1/2}
                    try:
                        self.dshpol = append(self.dshpol, [temp], axis=0)
                    except ValueError:
                        self.dshpol = [temp]
                    # Add the laser wavelength
                    self.b_e_frequencies = append(self.b_e_frequencies,
                                                  p.b_e_frequencies[0])
            else:
                #---------------------------------------------------------------
                # Zhongwei: In case of multiple groups of polarizability tensors
                temp = [p.polarizability[group] - m.polarizability[group]]
                #---------------------------------------------------------------
                # On first run self.polarizability is None, which raises
                # ValueError so we can initiallize the tensor before appending
                # to it.
                try:
                    self.qm_pol = append(self.qm_pol, temp,
                                                 axis=0)
                except ValueError:
                    self.qm_pol = temp
                #----------------------------------------------------------------------
                # Zhongwei: In case of multiple groups of polarizability tensors
                self.e_frequencies = append(self.e_frequencies, p.e_frequencies[group])
                #----------------------------------------------------------------------
            if gdipder and excitation is None:
                temp = (p.dipole - m.dipole)
                self.dgdip = append(self.dgdip, temp)
            # Remove the _LETTER for degenerate freq (if needed)
            tmp_num = num.rstrip('_abcdefghijklmnopqrstuvwxyz')
            vfreq = append(vfreq, float(tmp_num))

        # Reshape the ground state dipole derivatives, for convenience.
        if gdipder:
            self.dgdip = reshape(self.dgdip,(len(nums),3))
        # Reshape the transition dipole derivatives, for convenience.
        if tdipder:
            self.dtdip = reshape(self.dtdip,(len(nums),numexcite,3))
        # Reshape the transition quadrupole derivatives
        if len(self.TQM) != 0:
            self.TQM = reshape(self.TQM,(len(nums),numexcite,6))
        # Reshape the transition magnetic dipole derivatives
        if len(self.MDM) != 0 and 'CD SPECTRUM' in self.calctype:
            self.MDM = reshape(self.MDM,(len(nums),numexcite,3))

        # Eliminate modes that were not calculated.
        if self.nmodes > len(vfreq):
            from numpy import delete, array, where, zeros
            # Figure out the degeneracy of each frequency
            degeneracy = zeros((len(self.v_frequencies)), dtype=int)
            dn = 0
            lm = -1.0
            for i in range(len(self.v_frequencies)):
                if round(self.v_frequencies[i], 2) == round(lm, 2) or rnd(self.v_frequencies[i], 2) == rnd(lm, 2):
                    dn += 1
                else:
                    dn = 0
                degeneracy[i] = dn
                lm = self.v_frequencies[i]
            # Find whether that mode exists, otherwise delete it
            index = array([], dtype=int)
            sdegen = array(['','b','c','d','e','f','g','h','i'])
            for i in range(len(self.v_frequencies)):
                for j in range(len(vfreq)):
                    try:
                        sindex = where(sdegen==nums[j].split('_')[1])[0][0]
                    except IndexError:
                        sindex = 0
                    # Zhongwei: Test for the round-off issue
                    #if i == 0:
                    #   print(round(vfreq[j],2),rnd(vfreq[j],2))
                    #if j == 0:
                    #   print(round(self.v_frequencies[i],2),rnd(self.v_frequencies[i],2))
                    # Zhongwei: Fix the round-off issue
                    if (abs(round(self.v_frequencies[i], 2) - round(vfreq[j],2))
                        <= 0.01 and degeneracy[i] == sindex):
                    #if ((round(self.v_frequencies[i], 2) == round(vfreq[j], 2) or rnd(self.v_frequencies[i], 2) == rnd(vfreq[j], 2))
                    #    and degeneracy[i] == sindex):
                        break
                    else:
                        pass
                else:
                    index = append(index, i)
            self.IR = delete(self.IR, index)
            self.v_frequencies = delete(self.v_frequencies, index)
            self.normal_modes = delete(self.normal_modes, index, axis=0)
            self.nmodes = len(self.normal_modes)
            if excitation:
                self.es = excitation
            else:
                # Determine npol.
                self.npol = self.nmodes

        # Grab the stepsize for each normal mode
        if poop:
            sQ = ones(self.nmodes)*.01
        else:
            sQ = self.step_size(sR)
        if gdipder and excitation is None:
            for num in range(len(nums)):
                self.dgdip[num] = self.dgdip[num] / (2 * sQ[num])
        if excitation:
            from .constants import HBAR, LIGHT, PI, AMU, WAVENUM2HART, M2BOHR
            from .constants import BOHR2ANGSTROM as B2A
            from numpy import sqrt, array, einsum
            # Divide by the stepsize to make the derivatives
            self.deltas = self.deltas / ( 2 * sQ )
            # Multiply by a conversion to make dimentionless
            CONVERSION = sqrt(HBAR
                         / ( 2 * PI * LIGHT * 100 * self.v_frequencies * AMU ))
            self.deltas *= M2BOHR(CONVERSION)
            # Divide by frequency to make a dimentionless delta
            self.deltas = array(self.deltas / ( -self.v_frequencies * WAVENUM2HART ), dtype=float)

            # Add RAMAN to calctype
            self.calctype.update(['DELTAS', 'EXCITED STATE'])

            # Calculate the ground state dipole moment derivative
            if gdipder:
                for num in range(len(nums)):
                    self.dgdip[num] = self.dgdip[num] / (2 * sQ[num])
            # Calculate the transition dipole moment derivative
            if tdipder:
                self.dtdip = einsum('ij...,i...->ij...', self.dtdip, (1./sQ)) * 0.5
                self.TDM = self.dtdip # for consistency
                # Calculate the transition quadrupole moment derivatives
                if len(self.TQM) != 0: 
                    self.TQM = einsum('ij...,i...->ij...', self.TQM, (1./sQ)) * 0.5
                # Calculate the transition magnetic dipole derivatives
                if (len(self.MDM) != 0 and 'CD SPECTRUM' in self.calctype):
                    self.MDM = einsum('ij...,i...->ij...', self.MDM, (1./sQ)) * 0.5
        elif hpol:
            # Convert ds 1-d array to be an array of 3x3x3 tensors.
            from numpy import zeros
            temp = zeros((len(sQ),3,3,3))
            sQ = array([temp[i] + sQ[i] for i in range(len(sQ))])
            # Divide by the stepsize to make the derivatives
            self.dhpol = self.dhpol / ( 2 * sQ )
            # Add HYPER-RAMAN and HYPERPOLARIZABILITY to calctype
            self.calctype.update(['HYPER-RAMAN', 'HYPERPOLARIZABILITY'])
            if 'complex' in str(self.dhpol.dtype):
                self.calctype.add('FD')
            else:
                self.calctype.add('STATIC')
        elif shpol:
            # Convert ds 1-d array to be an array of 3x3x3x3 tensors.
            from numpy import zeros
            temp = zeros((len(sQ),3,3,3,3))
            sQ = array([temp[i] + sQ[i] for i in range(len(sQ))])
            # Divide by the stepsize to make the derivatives
            self.dshpol = self.dshpol / ( 2 * sQ )
            # Add SECOND HYPER-RAMAN and SECOND HYPERPOLARIZABILITY to calctype
            self.calctype.update(['SECOND HYPER-RAMAN', 'SECOND HYPERPOLARIZABILITY'])
            if 'complex' in str(self.dshpol.dtype):
                self.calctype.add('FD')
            else:
                self.calctype.add('STATIC')
        else:
            # Convert ds 1-d array to be an array of 3x3 tensors.
            # NumPy cannot broadcast a 1-d array onto a 3-d array.
            from numpy import zeros
            temp = zeros((len(sQ),3,3))
            sQ = array([temp[i] + sQ[i] for i in range(len(sQ))])
            # Divide by the stepsize to make the derivatives
            self.qm_pol = self.qm_pol / ( 2 * sQ )

            # Add RAMAN and polarizability to calctype
            self.calctype.update(['RAMAN', 'POLARIZABILITY'])
            if 'complex' in str(self.polarizability.dtype):
                self.calctype.add('FD')
            else:
                self.calctype.add('STATIC')


    def scale_vfreq(self, scale=1.0):
        '''Scales the vibrational frequency by a scaling factor. Default is unity.
        The following references is a good list of suggested scaling factors.

        J. Phys. Chem. A (2007), 111, 11683.

        '''
        self.v_frequencies *= scale

    def scattering_factor(self):
        '''Calculates the Raman scattering factor in units of
        :math:`\\frac{\AA{A}^4}{amu}`.

        Raises :py:exc:`AssertionError` if not a **RAMAN** calculation.

        '''
        from .constants import BOHR2ANGSTROM as B2A
        from numpy import absolute, zeros_like

        assert 'RAMAN' in self.calctype, ('scattering_factor(): '
                                                    'Not a RAMAN calculation.')

        # See if the scattering factors have already been collected
        # If not, calculate them now
        if self._raman is None:
            iso = self.isotropic()
            ani = self.anisotropic2()
            self._raman = absolute(
                             ( 45 * iso.conjugate() * iso + 7 * ani ) * B2A**4)

        # Now return the scattering factors
        return self._raman

    def raman_cross_section(self, laser=514.5, component='all', **kwargs):
        """Return the differential Raman cross section for each normal mode.

        The option *laser* is the laser wavelength to assume if the
        calculation was performed at the static limit.  The default is
        514.5.  The units are nm.

        To calculate the absolute differential raman scattering
        cross section:

        .. math::
           \\frac{d\\sigma}{d\\omega} = \\frac{h}{8c\\lambda_z\\epsilon_0^2}
           * \\frac{(\\lambda_0 - \\lambda_z)^4}
           {45(1-\\exp(\\frac{-hc\\lambda_z}{k_BT}))}*S

        where :math:`\\lambda_0` is the incident frequency in wavenumbers and
        :math:`\\lambda_z` is the vibrational frequency in wavenumbers and
        S is the :py:attr:`scattering_factor`.

        The squared polarizability derivatives should be in units of
        :math:`\\frac{C^2m^2}{V^2kg}`, but they are stored as squared
        derivatives of polarizability volumes in :math:`\\frac{\\AA^4}{amu}`.

        The conversion factor is:

        .. math:: \\frac{(4\\pi\\epsilon_0)^2*10^{-40}}{\\text{atomic mass unit}}

        Therefore, the total conversion factor is:

        .. math:: \\frac{2h\\pi^2*10^{40}}{\\text{atomic mass unit} * c}

        The cross section is returned in :math:`\\frac{cm^2}{sr}`.

        Raises :py:exc:`AssertionError` if not a **RAMAN** calculation.
        """
        from .constants import PI, PLANCK, LIGHT, BOLTZMAN, AMU
        from .constants import HART2WAVENUM, NM2WAVENUM, WAVENUM2INVM, M2CM
        from .constants import BOHR2ANGSTROM
        from numpy import exp, absolute

        assert 'RAMAN' in self.calctype, ('raman_cross_section(): '
                                                    'Not a RAMAN calculation.')

        CONVERSION = 2 * PI**2 * PLANCK * 1E-40 / ( LIGHT * AMU )
        EXPARG = PLANCK * LIGHT / BOLTZMAN
        TEMPERATURE = 298 # K
        # Default to given laser nm if calc was done at static limit
        if not self.e_frequencies[0]:
            print (laser)
            lambda_0 = NM2WAVENUM(laser)
        # Convert to wavenumber from hartrees if not at static limit.
        else:
            lambda_0 = HART2WAVENUM(self.e_frequencies[0])
        # Below, self.v_frequencies corresponds to lambda_z
        boltzfact = ( 1 - exp( - EXPARG * WAVENUM2INVM(self.v_frequencies) /
                                                                TEMPERATURE ) )
        freq = WAVENUM2INVM(lambda_0 - self.v_frequencies)**4
        # Convert from m^2 to cm^2

        # Use only one component to calculate raman cross sections -- Pengchong #
        if component=='zz':
            pol=self.polarizability[:,2,2]
            scat = absolute((45*pol.conjugate()*pol)*BOHR2ANGSTROM**4)
            print('Considering only ZZ component')
        elif component == 'xx':
            pol=self.polarizability[:,0,0]
            scat = absolute((45*pol.conjugate()*pol)*BOHR2ANGSTROM**4)
            print('Considering only XX component')
        elif component == 'yy': 
            pol=self.polarizability[:,1,1]
            scat = absolute((45*pol.conjugate()*pol)*BOHR2ANGSTROM**4)
            print('Considering only YY component')
        else:
            scat = self.scattering_factor()

        return M2CM(M2CM(( scat * CONVERSION ) * freq /
                       ( 45 * boltzfact * WAVENUM2INVM(self.v_frequencies) ) ))

    def hyperraman_cross_section(self, laser=1064, component='all', **kwargs):
        '''Return the differential hyper-Raman cross section for each normal mode.

        The option *laser* is the laser wavelength to assume if the
        calculation was performed at the static limit.  The default is
        1064 (fundamental of an Nd-YAG laser).  The units are nm.

        The differential hyper-Raman cross section is calculated using:

        .. math::
           \\frac{d\\sigma}{d\\Omega} = \\frac{2h^4\\alpha_f^3}{c^3e^6}
           * \\frac{(2\\lambda_0 - \\lambda_z)^4 * (2\\lambda'_0)}
           {\\lambda'_z(1-\\exp(\\frac{-hc\\lambda_z}{k_BT}))}
           *\\langle \\beta'_{ijk}\\rangle^2

        where :math:`\\lambda_0` is the incident frequency in Hertz,
        :math:`\\lambda_z` is the vibrational frequency in Hertz,
        :math:`\\lambda'_z` is the vibrational frequency in wavenumbers,
        and :math:`\\lambda'_0` is the incident frequency in Joules.

        The squared hyperpolarizability derivatives need to be in units of
        :math:`\\frac{C^6m^6}{J^4m^2kg}`, but they are stored as squared
        derivatives of hyperpolarizabilities
        :math:`\\frac{e^6a_0^6}{E_h^4\\AA^2amu}`.

        To convert the squared hyperpolarizability derivatives to SI units, the
        following conversion is used:

        .. math:: 1
                  \\left[\\frac{e^6a_0^6}
                               {E_h^4\\AA^2\\text{amu}}\\right]
                  =
                  2.210923386\\times 10^{-58}
                  \\left[\\frac{C^6m^6}
                               {J^4m^2kg}\\right]

        The cross section is returned in :math:`\\frac{cm^4 s}{photon sr}`.

        Raises :py:exc:`AssertionError` if not a **HYPER-RAMAN** calculation.
        '''

        from .constants import PLANCK, LIGHT, M2CM, BOLTZMAN, AMU
        from .constants import ELEM_CHARGE, HART2JOULE, BOHR2CM, ANGSTROM2M
        from .constants import FINE_STRUCTURE, CM2M, M2NM, WAVENUM2INVM
        from .constants import HART2WAVENUM
        from numpy import exp

        assert 'HYPER-RAMAN' in self.calctype, ('cross_section(): '
                                            'Not a HYPER-RAMAN calculation.')

        # The factor of 1E+8 in the numerator is used to convert the cross
        # section from m^4 s photon^{-1} sr^{-1} to cm^4 s photon^{-1} sr^{-1}.
        # I leave this in the conversion so that it is not confusing.

        CONVERSION1 = ELEM_CHARGE**6 * BOHR2CM**4 * 100**4
        CONVERSION1 = CONVERSION1/(100**4 * HART2JOULE**4 * AMU)

        CONVERSION2 = 2 * PLANCK**4 * FINE_STRUCTURE**3
        CONVERSION2 = CONVERSION2/(LIGHT**3 * M2CM * ELEM_CHARGE**6)

        # Exponential part of the Boltzmann distribution
        EXPARG = PLANCK * LIGHT / BOLTZMAN
        TEMPERATURE = 298 # K (to be consistent with the Raman)

        # Default to given laser nm if calc was done at static limit.
        if 'STATIC' in self.calctype:
            print (laser)
            lambda_0joule = 2*M2NM*CM2M/float(laser)              # cm^-1
            lambda_0hz = lambda_0joule*LIGHT/CM2M                 # s^-1
            lambda_0joule = lambda_0joule*HART2JOULE/HART2WAVENUM # J
        # Convert to wavenumber from hartrees if not at static limit.
        else:
            lambda_0hz = 2*HART2WAVENUM(self.b_e_frequencies)*LIGHT/CM2M
            lambda_0joule = 2*HART2JOULE(self.b_e_frequencies)

        # Below, self.v_frequencies corresponds to lambda_z
        boltzfact = ( 1 - exp( - EXPARG * WAVENUM2INVM(self.v_frequencies) /
                                                                TEMPERATURE ) )
        frq4th = (lambda_0hz - (self.v_frequencies)*LIGHT/CM2M)**4

        # Output the differential hyper-Raman cross section
        return ( CONVERSION1 * CONVERSION2 * frq4th * self.hpol_average(component=component) * lambda_0joule/
                 (boltzfact * self.v_frequencies))

    def secondhyperraman_cross_section(self, laser=1596, **kwargs):
        '''Return the differential second hyper-Raman cross section for each normal mode.

        The option *laser* is the laser wavelength to assume if the
        calculation was performed at the static limit.  The default is
        1596 (fundamental of an Nd-YAG laser).  The units are nm.

        The differential second hyper-Raman cross section is calculated using:

        .. math::
           \\frac{d\\sigma}{d\\Omega} = \\frac{8h^5\\alpha_f^4}{c^3e^8}
           * \\frac{(3\\lambda_0 - \\lambda_z)^4 * (3\\lambda'_0*\\lambda'_0)}
           {\\lambda'_z(1-\\exp(\\frac{-hc\\lambda_z}{k_BT}))}
           *\\langle \\gamma'_{ijkl}\\rangle^2

        where :math:`\\lambda_0` is the incident frequency in Hertz,
        :math:`\\lambda_z` is the vibrational frequency in Hertz,
        :math:`\\lambda'_z` is the vibrational frequency in wavenumbers,
        and :math:`\\lambda'_0` is the incident frequency in Joules.

        The squared hyperpolarizability derivatives need to be in units of
        :math:`\\frac{C^8m^8}{J^6m^2kg}`, but they are stored as squared
        derivatives of hyperpolarizabilities
        :math:`\\frac{e^8a_0^8}{E_h^6\\AA^2amu}`.

        To convert the squared second hyperpolarizability derivatives to SI units,
        the following conversion is used:

        .. math:: 1
                  \\left[\\frac{e^8a_0^8}
                               {E_h^6\\AA^2\\text{amu}}\\right]
                  =
                  8.361316595\\times 10^{-82}
                  \\left[\\frac{C^8m^8}
                               {J^6m^2kg}\\right]

        The cross section is returned in :math:`\\frac{cm^6 s^2}{photon^2 sr}`.

        Raises :py:exc:`AssertionError` if not a **SECOND HYPER-RAMAN** calculation.
        '''

        from .constants import PLANCK, LIGHT, M2CM, BOLTZMAN, AMU
        from .constants import ELEM_CHARGE, HART2JOULE, BOHR2CM, ANGSTROM2M
        from .constants import FINE_STRUCTURE, CM2M, M2NM, WAVENUM2INVM
        from .constants import HART2WAVENUM
        from numpy import exp

        assert 'SECOND HYPER-RAMAN' in self.calctype, ('cross_section(): '
                                            'Not a SECOND HYPER-RAMAN calculation.')

        # The factor of 1E+12 in the numerator is used to convert the cross
        # section from m^6 s^2 photon^{-2} sr^{-1} to cm^6 s^2 photon^{-2} sr^{-1}.
        # I leave this in the conversion so that it is not confusing.

        CONVERSION1 = ELEM_CHARGE**8 * BOHR2CM**6 * 100**6
        CONVERSION1 = CONVERSION1/(100**6 * HART2JOULE**6 * AMU)

        CONVERSION2 = 8 * PLANCK**5 * FINE_STRUCTURE**4
        CONVERSION2 = CONVERSION2/(LIGHT**3 * M2CM * ELEM_CHARGE**8)

        # Exponential part of the Boltzmann distribution
        EXPARG = PLANCK * LIGHT / BOLTZMAN
        TEMPERATURE = 298 # K (to be consistent with the Raman and hyper-Raman)

        # Default to given laser nm if calc was done at static limit.
        if 'STATIC' in self.calctype:
            print(laser)
            lambda_0joule = 3*M2NM*CM2M/laser                     # cm^-1
            lambda_0hz = lambda_0joule*LIGHT/CM2M                 # s^-1
            lambda_0joule = lambda_0joule*HART2JOULE/HART2WAVENUM # J
        # Convert to wavenumber from hartrees if not at static limit.
        else:
            lambda_0hz = 3*HART2WAVENUM(self.b_e_frequencies)*LIGHT/CM2M
            lambda_0joule = 3*HART2JOULE(self.b_e_frequencies)

        # Below, self.v_frequencies corresponds to lambda_z
        boltzfact = ( 1 - exp( - EXPARG * WAVENUM2INVM(self.v_frequencies) /
                                                                TEMPERATURE ) )
        frq4th = (lambda_0hz - self.v_frequencies*LIGHT/CM2M)**4

        # Output the differential second hyper-Raman cross section
        # The extra 3 in the numerator comes from lambda_0joule**2
        # See ''sechyperramancrs.f90" in the TDSPEC program
        return ( CONVERSION1 * CONVERSION2 * frq4th * self.shpol_average() * lambda_0joule**2/
                 (3 * boltzfact * self.v_frequencies))

    #def collect_tensor_derivatives(self, dir=None, poop=False, sR=0.01):
    # Zhongwei, provide the option to choose the 1st or 2nd or 3rd ... group of
    # polarizability tensors from higher-order property calculations, default
    # is the last group, i.e., "group = -1".
    def collect_tensor_derivatives(self, dir=None, poop=False, sR=0.01, group=-1, **kwargs):
        '''Collects the derivatives of the polarizability, A-tensor, optical rotation,
        C-tensor, G-tensor, beta tensor.
        Should be similar to the :py:attr:`raman_cross_section` function above.
        Note that these two routines should be merged in the future.'''

        from glob import glob
        from . import collect
        from natsort import natsort_key
        from numpy import append, array, where, ones, reshape, absolute, argmax, einsum,float64
        from numpy import round as rnd
        import os

        # Verify that we aren't doing something dumb.
        assert 'FREQUENCIES' in self.calctype, ('collect_raman_derivatives(): '
                         'File', self.filename, "is not of type 'FREQUENCIES'")

        # Grab the files in the directory of interest.
        if dir is None:
            dir = os.curdir
        pfiles = glob(os.path.join(dir, 'mode*-p.out'))
        mfiles = glob(os.path.join(dir, 'mode*-m.out'))

        # Make sure that the number of plus and minus files is the same
        assert len(mfiles) == len(pfiles), ('collect_raman_derivatives(): '
                    'Number of minus and plus direction files does not match.')

        # Sort the filenames 'naturally', i.e. in numerical order, not ASCII.
        pfiles.sort(key=natsort_key)
        mfiles.sort(key=natsort_key)

        # Make a list of just the numbers (as strings)
        nums = []
        for file in pfiles:
            name = os.path.split(file)[1]
            num = name[4:] # remove 'mode' from front
            nums.append(num.split('-')[0]) # Keep only the number part

        # Initialize some stuff
        vfreq = array([], dtype=float)
        self.e_frequencies  = array([], dtype=float)
        self.dgdip          = None
        self.quadrupole     = None
        self.qm_pol         = None
        self.dim_pol        = None
        self.ord            = None
        self.atensor        = None
        self.gtensor        = None
        self.deltas = array([], dtype=float)
        self.hyperpolarizability = None
        self.beta_ddq       = None
        self.beta_dqd       = None
        self.beta_qdd       = None
        self.beta_dqq       = None
        self.beta_qdq       = None
        self.beta_qqd       = None
        self.beta_qqq       = None
        # Zhongwei: Frequency for nonlinear properties
        self.b_e_frequencies  = array([], dtype=float)
        # Xing: hirshfeld polarizability
        self.hirsh_pol     = None
 
        # Loop over the numbers and collect the info into ChemData objects
        # Replace or augment some of the data.
        for num in nums:
            pname = 'mode' + '-'.join((num, 'p')) + '.out'
            mname = 'mode' + '-'.join((num, 'm')) + '.out'
            pfile = os.path.join(dir, pname)
            mfile = os.path.join(dir, mname)
            try:
                p = collect(pfile)
                m = collect(mfile)
            except IOError as e:
                sys.exit(str(e))

            # Find the difference between the tensors
            if p.polarizability is not None and m.polarizability is not None:
                #---------------------------------------------------------------
                # Zhongwei: In case of multiple groups of polarizability tensors
                temp = [p.polarizability[group] - m.polarizability[group]]
                #---------------------------------------------------------------
                # On first run self.polarizability is None, which raises
                # ValueError so we can initiallize the tensor before appending
                # to it.
                try:
                    self.qm_pol = append(self.qm_pol, temp, axis=0)
                except ValueError:
                    self.qm_pol = temp
                #----------------------------------------------------------------------
                # Zhongwei: In case of multiple groups of polarizability tensors
                self.e_frequencies = append(self.e_frequencies, p.e_frequencies[group])
                #----------------------------------------------------------------------
                self.subkey.add('POLARIZABILITY DERIVATIVES')
            # do dipole derivatives
            if p.dipole is not None and m.dipole is not None:
                temp = array(p.dipole-m.dipole)
                try:
                    self.dgdip = append(self.dgdip,temp,axis=0)
                except ValueError:
                    self.dgdip = array(temp)
                self.subkey.add("DIPOLE DERIVATIVES")
            # do quadrupole derivatives
            if p.quadrupole is not None and m.quadrupole is not None:
                temp = array(p.quadrupole - m.quadrupole)
                try:
                    self.quadrupole = append(self.quadrupole,temp,axis=0)
                except ValueError:
                    self.quadrupole = array(temp)
                self.subkey.add("QUADRUPOLE DERIVATIVES")
            # Do the same for the A-tensor
            if p.atensor is not None and m.atensor is not None:
                temp   = array(p.atensor - m.atensor)
                try:
                    self.atensor = append(self.atensor, temp, axis=0)
                except ValueError:
                    self.atensor = array(temp)
                self.subkey.add('A-TENSOR DERIVATIVES')
            # Do the same for the dip-mag dip G-tensor
            if p.gtensor is not None and m.gtensor is not None:
                temp = array(p.gtensor - m.gtensor)
                try:
                    self.gtensor = append(self.gtensor, temp, axis=0)
                except ValueError:
                    self.gtensor = array(temp)
                self.subkey.add('G-TENSOR DERIVATIVES')
            # Do the same for the quad-quad C-tensor
            if p.ctensor is not None and m.ctensor is not None:
                temp = array(p.ctensor - m.ctensor)
                try:
                    self.ctensor = append(self.ctensor, temp, axis=0)
                except ValueError:
                    self.ctensor = array(temp)
                self.subkey.add('C-TENSOR DERIVATIVES')
            # Do the same for the D-tensor
            if p.dtensor is not None and m.dtensor is not None:
                temp = array(p.dtensor - m.dtensor)
                try:
                    self.dtensor = append(self.dtensor, temp, axis=0)
                except ValueError:
                    self.dtensor = array(temp)
                self.subkey.add('D-TENSOR DERIVATIVES')
            # Do the same for the magnetizability
            if p.magnetizability is not None and m.magnetizability is not None:
                temp = array(p.magnetizability - m.magnetizability)
                try:
                    self.magnetizability = append(self.magnetizability, temp, axis=0)
                except ValueError:
                    self.magnetizability = array(temp)
                self.subkey.add('MAGNETIZABILITY DERIVATIVES')
            # Hyperpolarizability derivatives
            if p.hyperpolarizability is not None and m.hyperpolarizability is not None:
                try:
                    temp = array(p.hyperpolarizability - m.hyperpolarizability)
                except TypeError:
                    if 'STATIC' in p.hyperpolarizability.keys() and 'STATIC' in m.hyperpolarizability.keys():
                        temp = array([p.hyperpolarizability['STATIC'] - m.hyperpolarizability['STATIC']])
                        # Zhongwei: Add the static laser wavelength (0)
                        self.b_e_frequencies = append(self.b_e_frequencies, 0)
                    # Zhongwei: Added FD case here
                    elif 'FD' in p.hyperpolarizability.keys() and 'FD' in m.hyperpolarizability.keys():
                        temp = array([p.hyperpolarizability['FD'] - m.hyperpolarizability['FD']])
                        # Zhongwei: Add the laser wavelength
                        self.b_e_frequencies = append(self.b_e_frequencies, p.b_e_frequencies[0])
                try:
                    self.hyperpolarizability = append(self.hyperpolarizability, temp, axis=0)
                except ValueError:
                    self.hyperpolarizability = array(temp)
                self.subkey.add('HYPERPOLARIZABILITY DERIVATIVES')
            # Beta DDQ derivatives
            if p.beta_ddq is not None and m.beta_ddq is not None:
                temp = array(p.beta_ddq - m.beta_ddq)
                try:
                    self.beta_ddq = append(self.beta_ddq, temp, axis=0)
                except ValueError:
                    self.beta_ddq = array(temp)
                self.subkey.add('BETA_DDQ DERIVATIVES')
            # Beta DQD derivatives
            if p.beta_dqd is not None and m.beta_dqd is not None:
                temp = array(p.beta_dqd - m.beta_dqd)
                try:
                    self.beta_dqd = append(self.beta_dqd, temp, axis=0)
                except ValueError:
                    self.beta_dqd = array(temp)
                self.subkey.add('BETA_DQD DERIVATIVES')
            # Beta QDD derivatives
            if p.beta_qdd is not None and m.beta_qdd is not None:
                temp = array(p.beta_qdd - m.beta_qdd)
                try:
                    self.beta_qdd = append(self.beta_qdd, temp, axis=0)
                except ValueError:
                    self.beta_qdd = array(temp)
                self.subkey.add('BETA_QDD DERIVATIVES')
            # Beta DQQ derivatives
            if p.beta_dqq is not None and m.beta_dqq is not None:
                temp = array(p.beta_dqq - m.beta_dqq)
                try:
                    self.beta_dqq = append(self.beta_dqq, temp, axis=0)
                except ValueError:
                    self.beta_dqq = array(temp)
                self.subkey.add('BETA_DQQ DERIVATIVES')
            # Beta QDQ derivatives
            if p.beta_qdq is not None and m.beta_qdq is not None:
                temp = array(p.beta_qdq - m.beta_qdq)
                try:
                    self.beta_qdq = append(self.beta_qdq, temp, axis=0)
                except ValueError:
                    self.beta_qdq = array(temp)
                self.subkey.add('BETA_QDQ DERIVATIVES')
            # Beta QQD derivatives
            if p.beta_qqd is not None and m.beta_qqd is not None:
                temp = array(p.beta_qqd - m.beta_qqd)
                try:
                    self.beta_qqd = append(self.beta_qqd, temp, axis=0)
                except ValueError:
                    self.beta_qqd = array(temp)
                self.subkey.add('BETA_QQD DERIVATIVES')
            # Beta QQQ derivatives
            if p.beta_qqq is not None and m.beta_qqq is not None:
                temp = array(p.beta_qqq - m.beta_qqq)
                try:
                    self.beta_qqq = append(self.beta_qqq, temp, axis=0)
                except ValueError:
                    self.beta_qqq = array(temp)
                self.subkey.add('BETA_QQQ DERIVATIVES')
            # hirshfeld polarizability derivatives
            if p.hirsh_pol is not None and m.hirsh_pol is not None:
                temp = array(p.hirsh_pol[:self.natoms]-m.hirsh_pol[:self.natoms])
                try:
                    self.hirsh_pol = append(self.hirsh_pol,temp,axis=0)
                except ValueError:
                    self.hirsh_pol = array(temp)
                self.subkey.add('SHIRSHFELD_POLARIZABILITY DERIVATIVES')
            # Remove the _LETTER for degenerate freq (if needed)
            tmp_num = num.rstrip('_abcdefghijklmnopqrstuvwxyz')
            vfreq = append(vfreq, float(tmp_num))

        # Eliminate modes that were not calculated.
        if self.nmodes > len(vfreq):
            from numpy import delete, array, where, zeros
            # Figure out the degeneracy of each frequency
            degeneracy = zeros((len(self.v_frequencies)), dtype=int)
            dn = 0
            lm = -1.0
            for i in range(len(self.v_frequencies)):
                if (round(self.v_frequencies[i], 2) == round(lm, 2) or rnd(self.v_frequencies[i], 2) == rnd(lm, 2)):
                    dn += 1
                else:
                    dn = 0
                degeneracy[i] = dn
                lm = self.v_frequencies[i]
            # Find whether that mode exists, otherwise delete it
            index = array([], dtype=int)
            sdegen = array(['','b','c','d','e','f','g','h','i'])
            for i in range(len(self.v_frequencies)):
                for j in range(len(vfreq)):
                    try:
                        sindex = where(sdegen==nums[j].split('_')[1])[0][0]
                    except IndexError:
                        sindex = 0
                    #if ((round(self.v_frequencies[i], 2) == round(vfreq[j], 2) or rnd(self.v_frequencies[i], 2) == rnd(vfreq[j], 2))
                    #    and degeneracy[i] == sindex):
                    if (abs(round(self.v_frequencies[i], 2) - round(vfreq[j],2))
                        <= 0.01 and degeneracy[i] == sindex):
                        break
                    else:
                        pass
                else:
                    index = append(index, i)
            self.IR = delete(self.IR, index)
            self.v_frequencies = delete(self.v_frequencies, index)
            self.normal_modes = delete(self.normal_modes, index, axis=0)
            self.nmodes = len(self.normal_modes)
            # Determine npol.
            self.npol = self.nmodes

        # Grab the stepsize for each normal mode
        if poop:
            sQ = ones(self.nmodes)*.01
        else:
            sQ = self.step_size(sR)

        # Convert ds 1-d array to be an array of 3x3 tensors.
        # NumPy cannot broadcast a 1-d array onto a 3-d array.
        # Divide by the stepsize to make the derivatives
        fsQ = 1. / (2. * sQ)
        if 'DIPOLE DERIVATIVES' in self.subkey:
            self.dgdip = self.dgdip.reshape(len(self.e_frequencies),3)
            self.dgdip = array(self.dgdip,dtype=float64)
            self.dgdip = einsum('ia,i->ia', self.dgdip, fsQ)
        if 'QUADRUPOLE DERIVATIVES' in self.subkey:
            self.quadrupole = self.quadrupole.reshape(len(self.e_frequencies),6)
            self.quadrupole = array(self.quadrupole,dtype=float64)
            self.quadrupole = einsum('ia,i->ia', self.quadrupole, fsQ)
        if 'POLARIZABILITY DERIVATIVES' in self.subkey:
            self.qm_pol = einsum('iab,i->iab', self.qm_pol, fsQ)
        if 'G-TENSOR DERIVATIVES' in self.subkey:
            self.gtensor = einsum('iab,i->iab', self.gtensor, fsQ)
        if 'A-TENSOR DERIVATIVES' in self.subkey:
            self.atensor = einsum('iab,i->iab', self.atensor, fsQ)
        if 'C-TENSOR DERIVATIVES' in self.subkey:
            self.ctensor = einsum('iab,i->iab', self.ctensor, fsQ)
        if 'D-TENSOR DERIVATIVES' in self.subkey:
            self.dtensor = einsum('iab,i->iab', self.dtensor, fsQ)
        if 'MAGNETIZABILITY DERIVATIVES' in self.subkey:
            self.magnetizability = einsum('iab,i->iab', self.magnetizability, fsQ)
        if 'HYPERPOLARIZABILITY DERIVATIVES' in self.subkey:
            self.hyperpolarizability = einsum('iabc,i->iabc', self.hyperpolarizability, fsQ)
        if 'BETA_DDQ DERIVATIVES' in self.subkey:
            self.beta_ddq = einsum('iabc,i->iabc', self.beta_ddq, fsQ)
        if 'BETA_DQD DERIVATIVES' in self.subkey:
            self.beta_dqd = einsum('iabc,i->iabc', self.beta_dqd, fsQ)
        if 'BETA_QDD DERIVATIVES' in self.subkey:
            self.beta_qdd = einsum('iabc,i->iabc', self.beta_qdd, fsQ)
        if 'BETA_DQQ DERIVATIVES' in self.subkey:
            self.beta_dqq = einsum('iabc,i->iabc', self.beta_dqq, fsQ)
        if 'BETA_QDQ DERIVATIVES' in self.subkey:
            self.beta_qdq = einsum('iabc,i->iabc', self.beta_qdq, fsQ)
        if 'BETA_QQD DERIVATIVES' in self.subkey:
            self.beta_qqd = einsum('iabc,i->iabc', self.beta_qqd, fsQ)
        if 'BETA_QQQ DERIVATIVES' in self.subkey:
            self.beta_qqq = einsum('iabc,i->iabc', self.beta_qqq, fsQ)
        if 'SHIRSHFELD_POLARIZABILITY DERIVATIVES' in self.subkey:
            self.hirsh_pol = self.hirsh_pol.reshape(len(self.e_frequencies),self.natoms,3,3)
            self.hirsh_pol = einsum('iabc,i->iabc', self.hirsh_pol, fsQ)

        # Add VROA to calctype
        if 'OPTICAL ROTATION' in self.subkey and 'A-TENSOR' in self.subkey:
            self.calctype.update(['VROA', 'POLARIZABILITY', 'A-TENSOR',
                                  'OPTICAL ROTATION', 'G-TENSOR'])
        elif 'OPTICAL ROTATION' in self.subkey:
            self.calctype.update(['POLARIZABILITY', 'OPTICAL ROTATION', 'G-TENSOR'])
        elif 'A-TENSOR' in self.subkey:
            self.calctype.update(['POLARIZABILITY', 'A-TENSOR'])
        if 'complex' in str(self.polarizability.dtype):
            self.calctype.add('FD')
        else:
            self.calctype.add('STATIC')

    def vroa_cross_section(self, laser=514.5):
        '''See raman_cross_section above'''
        from .constants import PI, PLANCK, LIGHT, BOLTZMAN, AMU
        from .constants import HART2WAVENUM, NM2WAVENUM, WAVENUM2INVM, M2CM
        from numpy import exp

        CONVERSION = PI**2 * PLANCK * 1E-40 / ( LIGHT * AMU )
        EXPARG = PLANCK * LIGHT / BOLTZMAN
        TEMPERATURE = 298 # K
        # Default to given laser nm if calc was done at static limit
        if not self.e_frequencies[0]:
            laser=raw_input("Please give the laser wavelength [514.5]: ") or 514.5
            lambda_0 = NM2WAVENUM(laser)
        # Convert to wavenumber from hartrees if not at static limit.
        else:
            lambda_0 = HART2WAVENUM(self.e_frequencies[0])
        # Below, self.v_frequencies corresponds to lambda_z
        boltzfact = ( 1 - exp( - EXPARG * WAVENUM2INVM(self.v_frequencies) /
                                                                TEMPERATURE ) )
        freq = WAVENUM2INVM(lambda_0 - self.v_frequencies)**4
        # Convert from m^2 to cm^2
        self.vroa_intensities['180deg'] = M2CM(M2CM(( self.vroa_intensities['180deg'] * CONVERSION ) * freq /
                           ( 45 * boltzfact * WAVENUM2INVM(self.v_frequencies) ) ))


    def calc_roa_intensities(self,  angle='180deg'):
        '''Calculates the ROA intensities based on the equation in ADF output:
        :math:`I^R(180^o)-I^L(180^o) = 48/c * [B(G)^2 + (B(A)^2)/3]`
        for the :math:`180^o` collection angle (back scattered).'''

        import numpy as np
        from numpy import zeros, array

        # Temporarily store the tensor derivatives
        alpha   = np.copy(self.polarizability)
        atensor = np.copy(self.atensor)
        optrot  = np.copy(self.ord)

        # Import the Levi-Civita tensor
        from .constants import LEVICIVITA3 as epsilon

        # figure out if complex tensors are needed
        lifetime = isinstance(alpha[0][0][0], complex)

        # Calculate B(G)**2
        BG  = zeros( (self.nmodes), dtype=complex)
        for i in range(self.nmodes):
            for a in range(3):
                for b in range(3):
                    BG[i] += (0.5*(3.0*alpha[i][a][b]*optrot[i][a][b].conjugate()
                             - alpha[i][a][a]*optrot[i][b][b].conjugate()))*1j
        # for some reason, the values above are 2000x larger:
        BG = BG / 2.0

        # Calculate B(A)**2
        BA  = zeros( (self.nmodes), dtype=complex)
        for i in range(self.nmodes):
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        for d in range(3):
                            if d == 0:
                                e = b
                            elif b == 0:
                                e = d
                            else:
                                e = b + d + 1
                            BA[i] += 0.5*(alpha[i][a][b]
                                        *epsilon[a][c][d]*atensor[i][c][e].conjugate()
                                        *self.e_frequencies[0])
        # again, values above are 2000x larger:
        BA = BA / 2.0

        # Calculate (aG')**2
        aG = zeros( (self.nmodes), dtype=complex)
        for i in range(self.nmodes):
            for a in range(3):
                for b in range(3):
                    aG[i] += (1.0/9.0)*1j*(alpha[i][a][a]
                             * optrot[i][b][b].conjugate())
        # again, values above are 2x larger
        aG = aG / 2.0

        # Print B(G)**2 and B(A)**2 if needed
        if True:
            if lifetime: print ('Real Part:')
            print (' Freq         (aG)**2       B(G)**2       B(A)**2   (x10E6 Ang**4/amu)')
            for i in range(self.nmodes):
                print ('{0:7.2f} {1:13.4f} {2:13.4f} {3:13.4f}'.format(self.v_frequencies[i],
                       aG[i].imag*1e3, BG[i].imag*1e3, BA[i].real*1e3))
            print ('')
            if lifetime:
                print ('Imaginary Part:')
                print (' Freq         (aG)**2       B(G)**2       B(A)**2   (x10E6 Ang**4/amu)')
                for i in range(self.nmodes):
                    print ('{0:7.2f} {1:13.4f} {2:13.4f} {3:13.4f}'.format(self.v_frequencies[i],
                           aG[i].real*1e3, BG[i].real*1e3, BA[i].imag*1e3))
                print ('')

        # Calculate Intensity for backscattered ROA
        if angle != '180deg': print ('WARNING: calculating for backscattered '
                                     'ROA anyway!')
        intensity = zeros ( (self.nmodes), dtype=float)
        for i in range(self.nmodes):
            intensity[i] = 48.0 * ( BG[i].imag + (BA[i].real)/3.0 )

        # Add intensity to appropriate variable
        self.vroa_intensities = {'180deg': intensity}
        if True:
        # multiply by kp
            self.vroa_cross_section(self)
        else:
            self.vroa_intensities[angle] = intensity


    def integrateSpec(self, property='cross section', md1=1, md2=None):
        '''Return the integration of a vibrational spectroscopic property.
       *md1* is the lowest mode to include in the sum

        *md2* is the highest mode to include in the sum, defaulting to
        py:attr:`nmodes`.

        *property* currently has three choices implemented:

         - 'cross section'
         - 'scattering factor'
         - 'IR'

        '''
        if property=='IR': assert 'FREQUENCIES' in self.calctype, (
                                'integrate(): Not a FREQUENCIES calculation.')

        # Convert mode number to index number
        md1 -= 1
        # Default to all.
        if md2 is None: md2 = self.nmodes
        if property == 'cross section':
            return self.cross_section()[md1:md2].sum()
        elif property == 'scatterting factor':
            return self.scattering_factor()[md1:md2].sum()
        elif property == 'IR':
            return self.IR[md1:md2].sum()
        else:
            raise ChemDataError ('integrate(): Invalid property: '+property)


