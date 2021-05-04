#! /usr/bin/env python
from __future__ import print_function

def raman_cross_section(x, y, e_frequencies=None, laser=514.5, **kwargs):
    '''This function is identical to the raman_cross_section function
    with the exception that the scattering factors are given by "y"
    instead of calculated from the polarizabilities.'''

    from ..constants import PI, PLANCK, LIGHT, BOLTZMAN, AMU
    from ..constants import HART2WAVENUM, NM2WAVENUM, WAVENUM2INVM, M2CM
    from numpy import exp

    CONVERSION = 2 * PI**2 * PLANCK * 1E-40 / ( LIGHT * AMU )
    EXPARG = PLANCK * LIGHT / BOLTZMAN
    TEMPERATURE = 298 # K
    # Zhongwei: Make sure the static case is working appropriately
    # Default to given laser nm if calc was done at static limit
    if not e_frequencies[0]:
        print (laser)
        lambda_0 = NM2WAVENUM(laser)
    # Convert to wavenumber from hartrees if not at static limit.
    else:
        lambda_0 = HART2WAVENUM(e_frequencies[0])
    # Below, self.v_frequencies corresponds to lambda_z
    boltzfact = ( 1 - exp( - EXPARG * WAVENUM2INVM(x) / TEMPERATURE ) )
    freq = WAVENUM2INVM(lambda_0 - x)**4
    # Convert from m^2 to cm^2
    return M2CM(M2CM(( y * CONVERSION ) * freq / ( boltzfact * WAVENUM2INVM(x) ) ))

def cars_cross_section(x, y, e_frequencies=None, laser=514.5, **kwargs):
    '''This function is identical to the raman_cross_section function
    with the exception that the scattering factors are given by "y"
    instead of calculated from the polarizabilities.'''

    from ..constants import ELEM_CHARGE, BOHR2ANGSTROM, HART2JOULE, AMU
    from ..constants import PI, FINE_STRUCTURE, PLANCK, LIGHT, BOLTZMAN
    from ..constants import VACUUM_EL, NM2HART, WAVENUM2HART
    from ..constants import HART2WAVENUM, NM2WAVENUM, WAVENUM2INVM, M2CM
    from numpy import exp

    # Zhongwei: Make sure the static case is working appropriately
    # Default to given laser nm if calc was done at static limit
    if not e_frequencies[0]:
        print (laser)
        lambda_0 = NM2WAVENUM(laser)
        omega_0 = HART2JOULE(NM2HART(lambda_0))
    # Convert to wavenumber from hartrees if not at static limit.
    else:
        lambda_0 = HART2WAVENUM(e_frequencies[0])
        omega_0 = HART2JOULE(e_frequencies[0])

#    CONVERSION1 = (( 4 * PI * VACUUM_EL )**2 * BOHR2ANGSTROM**4 * 1E-40 / ( AMU * HART2JOULE ))**2
    CONVERSION1 = ELEM_CHARGE**8 * BOHR2ANGSTROM**4 * 1E-40 / ( HART2JOULE**6 * AMU**2 )
    CONVERSION2 = 64 * PI**2 * FINE_STRUCTURE**4 * PLANCK**4 * LIGHT**2 / ELEM_CHARGE**8
    CONVERSION3 = ( PLANCK / ( 8 * PI**2 * LIGHT ))**2
    CONVERSION4 = omega_0 * ( omega_0 - HART2JOULE(WAVENUM2HART(x)) )
    CONVERSION = CONVERSION1 * CONVERSION2 * CONVERSION3 * CONVERSION4
    EXPARG = PLANCK * LIGHT / BOLTZMAN
    TEMPERATURE = 298 # K
    # Below, self.v_frequencies corresponds to lambda_z
    boltzfact = ( 1 - exp( - EXPARG * WAVENUM2INVM(x) / TEMPERATURE ) )
    frq4th = WAVENUM2INVM(lambda_0 + x)**4
    # Convert from m^6 to cm^6
    return 1E12 * (( y * CONVERSION ) * frq4th / ( boltzfact * WAVENUM2INVM(x)**2 ) )


def hyperraman_cross_section(vfreq, hpol, b_e_frequencies=None, laser=1064, component='all', **kwargs):
    '''Function identical to the chem package function of the same name.'''

    from ..constants import PLANCK, LIGHT, M2CM, BOLTZMAN, AMU
    from ..constants import ELEM_CHARGE, HART2JOULE, BOHR2CM, ANGSTROM2M
    from ..constants import FINE_STRUCTURE, CM2M, M2NM, WAVENUM2INVM
    from ..constants import HART2WAVENUM
    from numpy import exp, sqrt

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

    # Zhongwei: Make sure the static case is working appropriately
    #           Take the frequency-dependent case into account
    # Default to given laser nm if calc was done at static limit.
    if not b_e_frequencies[0]:
        print (laser)
        lambda_0joule = 2*M2NM*CM2M/float(laser)              # cm^-1
        lambda_0hz = lambda_0joule*LIGHT/CM2M                 # s^-1
        lambda_0joule = lambda_0joule*HART2JOULE/HART2WAVENUM # J
    # Convert to wavenumber from hartrees if not at static limit.
    else:
        lambda_0hz = 2*HART2WAVENUM(b_e_frequencies[0])*LIGHT/CM2M
        lambda_0joule = 2*HART2JOULE(b_e_frequencies[0])

    # Below, vfreq corresponds to lambda_z
    boltzfact = ( 1 - exp( - EXPARG * WAVENUM2INVM(vfreq) /
                                                            TEMPERATURE ) )
    frq4th = (lambda_0hz - vfreq*LIGHT/CM2M)**4

    # Zhongwei: in case we only want certain components for hyper-Raman/SEHRS
    if component != 'all':
       if len(component) == 3:
          print('Considering only ' + component.upper() + ' component')
          sym_change = { 'X' : '0', 'Y' : '1', 'Z' : '2', }
          sym_index = []
          for i in range(len(component)):
              sym_index.append(int(sym_change[component[i].upper()]))
          for i in range(len(hpol)):
              for j in range(3):
                  for k in range(3):
                      for l in range(3):
                          if j == sym_index[0] and k == sym_index[1] and l == sym_index[2]:
                             pass
                          else:
                             hpol[i,j,k,l] = 0.0
       else:
          print('The component ' + component.upper() + ' is invalid for hyperpolarizability.')
          print('Please use the XYZ format.')
          import sys
          sys.exit()

    # Output the differential hyper-Raman cross section
    hpol_avg = hpol_average(hpol)

    ##########################################################################
    # Zhongwei: print out the \beta components values for analyzing the
    # orientation information from SEHRS
    #for i in range(len(hpol)):
    #    print ('Mode:', vfreq[i])
    #    for j in range(3):
    #        for k in range(3):
    #            for l in range(3):
    #                print (j,k,l,hpol[i][j][k][l].real,hpol[i][j][k][l].imag)
    #    print ('-----------------------------------------------------------')
    # The ZZZ component only!
    #for i in range(len(hpol)):
    #    print (vfreq[i], sqrt((hpol[i][2][2][2].real)**2+(hpol[i][2][2][2].imag)**2))
    #print ('-----------------------------------------------------------')
    #########################################################################

    return ( CONVERSION1 * CONVERSION2 * frq4th * hpol_avg * lambda_0joule/
             (boltzfact * vfreq))

def hpol_average(tn):
    '''Copied from the chem package function of the same name.'''

    from numpy import zeros

    # Defining the coordinates
    first = range(3)
    second = range(3)
    third = range(3)

    cj = tn.conjugate()

    # Arrays storing the total hyperpolarizability averages.
    
    bmean = zeros((len(tn),2),dtype=complex)
    hpol_avg = zeros(len(tn), dtype=complex)
    b_avg = zeros((len(tn),15), dtype=complex)
    for n in xrange(len(tn)):
        b_avg[n][0] = sum([tn[n][i][i][i]*cj[n][i][i][i]
                               for i in first])

        b_avg[n][1] = sum([tn[n][i][i][j]*cj[n][i][i][j]
                               if i != j else 0.0 for i in first for j in second])

        b_avg[n][2] = sum([tn[n][i][i][i]*cj[n][i][j][j]
                               if i != j else 0.0 for i in first for j in second])

        b_avg[n][3] = sum([tn[n][j][i][i]*cj[n][i][i][j]
                               if i != j else 0.0 for i in first for j in second])

        b_avg[n][4] = sum([tn[n][i][i][i]*cj[n][j][j][i]
                               if i != j else 0.0 for i in first for j in second])

        b_avg[n][5] = sum([tn[n][j][i][i]*cj[n][j][i][i]
                               if i != j else 0.0 for i in first for j in second])

        b_avg[n][6] = sum([tn[n][i][i][j]*cj[n][j][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

        b_avg[n][7] = sum([tn[n][j][i][i]*cj[n][j][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

        b_avg[n][8] = sum([tn[n][i][i][j]*cj[n][k][k][j]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

        b_avg[n][9] = sum([tn[n][i][j][k]*cj[n][i][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

        b_avg[n][10] = sum([tn[n][i][j][k]*cj[n][j][i][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

        b_avg[n][11] = sum([tn[n][i][j][j]*cj[n][i][j][j]
                               if i != j else 0.0
                               for i in first for j in second])

        b_avg[n][12] = sum([tn[n][i][i][j]*cj[n][j][i][i]
                               if i != j else 0.0
                               for i in first for j in second])

        b_avg[n][13] = sum([tn[n][i][j][j]*cj[n][i][k][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

        b_avg[n][14] = sum([tn[n][i][i][k]*cj[n][j][j][k]
                               if (i != j and i != k and j != k) else 0.0
                               for i in first for j in second for k in third])

        bmean[n][0] = ( b_avg[n][0]/7 + 4*b_avg[n][1]/35 + 2*b_avg[n][2]/35
                      + 4*b_avg[n][3]/35 + 4*b_avg[n][4]/35 + b_avg[n][5]/35
                      + 4*b_avg[n][6]/105 + b_avg[n][7]/105 + 4*b_avg[n][8]/105
                      + 2*b_avg[n][9]/105 + 4*b_avg[n][10]/105)
        bmean[n][1] = ( b_avg[n][0]/35 + 4*b_avg[n][2]/105 - 4*b_avg[n][4]/70
                      + 8*b_avg[n][1]/105 + 3*b_avg[n][11]/35 - 4*b_avg[n][12]/70
                      + b_avg[n][13]/35 - 4*b_avg[n][14]/210 - 4*b_avg[n][6]/210
                      + 2*b_avg[n][9]/35 - 4*b_avg[n][10]/210)

        hpol_avg[n] = bmean[n][0] + bmean[n][1]

    return hpol_avg.real

def plot(Freq, Intensity, fwhm=20, minx=500, maxx=1700, step=1, sticks=True,
    fs=20, lw=1.5, noplot=False, spectrum=None, **kwargs):
    '''Generates a lorentzian from input Freq and Intensity and plots the resulting
    spectrum.'''

    from matplotlib import pyplot
    from .functions import sum_lorentzian
    from numpy import arange
    from math import pi, log10

    ckwargs = check_kwargs
    rkwargs = return_kwargs

    # Set some defaults
    axfs = int(fs*0.8)
    x = arange(minx,maxx,step)

    # Return if no plotting with sticks
    if noplot and sticks: return Freq, Intensity

    # Generate lorentzian
    if spectrum == 'CARS':
        x = Freq[:]
        y = Intensity[:]
    else:
        y = sum_lorentzian(x, Freq, Intensity, fwhm=fwhm)

    # Return arrays if no plotting
    if noplot: return x, y

    # Add plot to figure
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    scale = rkwargs('scale', **kwargs)
    if not scale:
        if spectrum == 'CID':
            if abs(y.min()) < 1e-3 and abs(y.max()) < 1e-3:
                scale = 1e3
            else:
                scale = 1.
        elif spectrum == 'CARS':
            if abs(y.max()) < 1e-92:
                cale = 1e95
            elif abs(y.max()) < 1e-87:
                scale = 1e90
            else:
                scale = 1e85
        elif spectrum == 'HRS' or spectrum == 'HYPERRAMAN' or spectrum == 'HSERS':
            if abs(y.max()) < 1e-62:
                scale = 1e65
            elif abs(y.max()) < 1e-59:
                scale = 1e62
            else:
                scale = 1e60
        elif abs(y.max()) < 1e-34 and abs(y.min()) < 1e-34:
            scale = 1e34
        else:
            scale = 1e30
    if sticks and not spectrum == 'CARS':
        for i in range(len(Freq)):
            ax.plot((Freq[i], Freq[i]), (0, Intensity[i]*2*(scale)/(fwhm*pi)), 'r', lw=lw)
    ax.plot(x,y*scale, ls='-', lw=lw, color='k')

    pyplot.xlim((minx,maxx))
    pyplot.rcParams['font.size'] = fs
    pyplot.xlabel('wavenumber ($cm^{-1}$)', fontsize=fs)

    # Figure out whether we're plotting CID or cross-section
    if spectrum == 'CID':
        if ckwargs('geom', 'FORWARD', **kwargs):
            string = '0^o'
        elif ckwargs('geom', 'NINETY_X', **kwargs):
            string = '90^o_x'
        elif ckwargs('geom', 'NINETY_Z', **kwargs):
            string = '90^o_z'
        else:
            string = '180^o'
        if abs(y.min()) < 1e-3 and abs(y.max()) < 1e-3:
            string2 = '$\\left( \\times 10^{-3} \\right)$'
        else:
            string2 = ''
        pyplot.ylabel('CID $\\left('+string+'\\right)$ '+string2, fontsize=fs)
        # Make sure MIN range is between -1e-5 and +1e-5
        if abs(y.min())*scale < 1e-3 and abs(y.max())*scale < 1e-3:
            pyplot.ylim((-1e-3,1e-3))
    elif spectrum == 'CARS':
        string = '{0}'.format(int(log10(1./scale)))
        pyplot.ylabel('$d\sigma/d\Omega \\left( 10^{'+string+'} '
                      ' \\frac{cm^6 s^2}{{photon}^2 sr}\\right)$', fontsize=fs)
    elif spectrum == 'HRS' or spectrum == 'HYPERRAMAN' or spectrum == 'HSERS':
        string = "{0}".format(int(log10(1./scale)))
        pyplot.ylabel("$d\sigma/d\Omega$ / $\\times 10^{"+string+"} \\frac{cm^4 s}{photon sr}$", fontsize=fs)
    else:
        string = '{0}'.format(int(log10(1./scale)))
        pyplot.ylabel('$d\sigma/d\Omega$ $\\left(10^{'+string+'} cm^2/sr\\right)$', fontsize=fs)
        # Make sure that range is between -1e-30 and +1e-30
        if abs(y.max())*scale < 1e-5 and y.min() >= 0.:
            pyplot.ylim((0.,1e-5))
        elif abs(y.max())*scale < 1e-5 and abs(y.min())*scale < 1e-5:
            pyplot.ylim((-1e-5,1e-5))

    pyplot.xticks(fontsize=axfs)
    pyplot.yticks(fontsize=axfs)

    # Are inverted axes needed?
    #if rkwargs('invertx', True, **kwargs):
    #    pyplot.gca().invert_xaxis()
    if rkwargs('invertx', False, **kwargs):
        pyplot.gca().invert_xaxis()
    if rkwargs('inverty', False, **kwargs):
        pyplot.gca().invert_yaxis()

    # Save figure if requested
    if ckwargs('savefig', **kwargs):
        filename = rkwargs('filename', 'default', **kwargs)
        filetype = rkwargs('format', 'png', **kwargs)
        if len(filename) < 5 or filename[-4] != '.':
            filename = filename + '.' + filetype
        dpi = rkwargs('dpi', 300., **kwargs)
        pyplot.savefig(filename, dpi=dpi, format=filetype, transparent=True)
        print ('Figure {0} saved!'.format(filename))

    # Show figure
    pyplot.show()


def print_field(E, FG=None):
    '''Pretty prints the E-field and gradient.'''

    Exyz = ['x','y','z']
    FGxyz = [['xx','xy','xz'],['yx','yy','yz'],['zx','zy','zz']]
    print ('')
    print ('Field and field gradient used are:')
    for i in range(3):
        for j in range(3):
            print ('E('+Exyz[i]+')_'+Exyz[j]+'   {0: 10.4f} {1: 10.4f}j'.format(E[i][j].real, E[i][j].imag))
    if FG is not None:
        for i in range(3):
            for j in range(3):
                for k in range(j,3):
                    print('FG('+Exyz[i]+')_'+FGxyz[j][k]+' {0: 10.4f} {1: 10.4f}j'.format(FG[i][j][k].real,
                                                                      FG[i][j][k].imag))
    print ('')


def return_dim_field(d, com, tr=None, debug=True, index=0,screenfac=1.0, polarization='ALL', gradient=True):
    '''Returns the local field and gradient from a DIM particle at 
    a given the center of mass.'''

    from ..  import collect
    from ..constants import ANGSTROM2BOHR, KRONECKER3 as dt
    from numpy import array, zeros, einsum, sqrt, pi, exp
    from copy import deepcopy
    from math import erf

    # calculate the local field at the molecule's center of mass
    if tr is None: tr = zeros((3)) # incase we need to translate the DIM system
    try:
        d = collect(d)
    except AttributeError:
        pass
    except TypeError:
        pass
    except IOError:
        if debug: print ('DIM system not present, using zero field and gradients...')
        return zeros((3,3)), zeros((3,3,3))
    dip = deepcopy(d.dim_dipoles['FD scattered'][index][:])
    dip = einsum('ijk->jik', dip)

    # account for the requested incident polarization
    pol = array([1.,1.,1.])
    if not ('X' in polarization.upper() or polarization.upper() == 'ALL'):
        pol[0] = 0.
    if not ('Y' in polarization.upper() or polarization.upper() == 'ALL'):
        pol[1] = 0.
    if not ('Z' in polarization.upper() or polarization.upper() == 'ALL'):
        pol[2] = 0.
    dip = einsum('ijk,j->ijk', dip, pol)

    R   = deepcopy(d.dim_coordinates[:])
    R   = array(R) * ANGSTROM2BOHR
    r   = array(com - R - tr, dtype=float)
    r1  = sqrt((r**2).sum(axis=1))
    r2  = r1*r1
    r3  = r2*r1
    r5  = r3*r2
    if gradient: r7  = r5*r2

    if screenfac == 1.0:
        screen = 3.3208 # this value is sqrt(1.447**2 + 1**2) converted to bohr (default)
    elif screenfac is None:
        screen = 1. # No screening
    else:
        screen = sqrt(7.4529 + ANGSTROM2BOHR(screenfac)**2) # use screenfac probe distance
    inr = 1. / screen
    s   = r1 * inr
    sf0 = zeros((len(s)), dtype=float)
    for i in range(len(sf0)):
        sf0[i] = erf(s[i])
    sf1 = ( 2. * s / sqrt(pi) ) * exp(-s**2)
    sf2 = ( 4. * inr**3 / sqrt(pi) ) * exp(-s**2)
    sf  = sf0 - sf1
    if gradient:
        sf3 = ( 4. * s**2 / ( sqrt(pi) * r2 )) * exp(-s**2)
        sf4 = ( 4. * s**3 / sqrt(pi) ) * exp(-s**2)

    E = 3. * einsum('i,ib,ic,iac,i->ab', sf, r, r, dip, (1/r5))
    E -= einsum('i,bc,iac,i', sf, dt, dip, (1/r3))
    E -= einsum('i,ib,ic,iac,i', sf2, r, r, dip, (1/r2))

    if gradient:
        FG  = -15 * einsum('i,ib,ic,id,iad,i->abc', sf, r, r, r, dip, (1/r7))
        FG += 3 * einsum('i,bc,id,iad,i->abc', sf, dt, r, dip, (1/r5))
        FG += 3 * einsum('i,bd,ic,iad,i->abc', sf, dt, r, dip, (1/r5))
        FG += 3 * einsum('i,cd,ib,iad,i->abc', sf, dt, r, dip, (1/r5))

        FG += 3 * einsum('i,ib,ic,id,iad,i->abc', sf3, r, r, r, dip, (1/r5))
        FG -= einsum('i,bc,id,iad,i->abc', sf3, dt, r, dip, (1/r3))

        FG += 3 * einsum('i,ib,id,ic,iad,i->abc', sf3, r, r, r, dip, (1/r5))
        FG -= einsum('i,bd,ic,iad,i->abc', sf3, dt, r, dip, (1/r3))

        FG -= einsum('i,ib,cd,iad,i->abc', sf4, r, dt, dip, (1/r5))
        FG -= einsum('i,ib,ic,id,iad,i->abc', sf4, r, r, r, dip, (1/r7))
        FG += 2 * einsum('i,i,ib,ic,id,iad,i->abc', sf4, s**2, r, r, r, dip, (1/r7))

    # return fields
    if gradient:
        return E, FG
    else:
        return E


def generate_field (f, a=None, r=None, alpha=None, ei=None, SFG=False, freq_VIS=354):
    '''Generates the field and field gradient at the center of mass of
    the molecule (cm) given the the vector from the NP to the molecule (r)
    and the polarizability of the nanoparticle (alpha).
    e is the pertubating field vector, and is taken as 1 if no fields given.'''
    
    from ..constants import NM2WAVENUM, WAVENUM2HART, NM2HART, ANGSTROM2BOHR 
    from ..constants import KRONECKER3 as dt
    from numpy import array, zeros, einsum
    from math import sqrt

    # modify r so that it describes the separation of your origins
    
    # calculate r and a3
    if r is None: r = array([35.,0.,0.])
    r = ANGSTROM2BOHR(r) #convert to bohr 
    #Zhongwei's note: The change above is made based on Xing's modification.
    #                 For DIM nanoparticle, this should be already in Bohr.
    #                 For isotropic sphere, this might not matter.
    r1 = sqrt( ( r * r.conj() ).sum().real )
    r3 = r1 * r1 * r1
    r5 = r3 * r1 * r1
    r7 = r5 * r1 * r1
    r9 = r7 * r1 * r1
    if a is None: a = r1 # a is the radius of the sphere
    else: a = a * ANGSTROM2BOHR # See Zhongwei's note above
    a3 = a**3

    if SFG:
        di_freqIR = zeros((len(f.v_frequencies)), dtype=float)
        di_freqSCAT = zeros((len(f.v_frequencies)), dtype=float)
        di_freqVIS = array([NM2WAVENUM(freq_VIS),NM2WAVENUM(550)]) #2nd term is a dummy so it's iterable
        di_path = '/usr/repo/applications/dim/dim/dielectrics/Ag_jc'
        for x in range(0, len(f.v_frequencies)):
            di_freqIR[x]  = f.v_frequencies[x] 
            di_freqSCAT[x] = f.v_frequencies[x] + di_freqVIS[0]
        # Change to Hartrees for spline
        di_freqSCAT = WAVENUM2HART(di_freqSCAT)
        di_freqVIS = WAVENUM2HART(di_freqVIS)
        di_freqIR  = WAVENUM2HART(di_freqIR)
    
        # Pull the dielectric function
        ei_VIS = expdie(di_path, di_freqVIS)
        ei_SCAT = expdie(di_path, di_freqSCAT)
        ei_IR  = expdie(di_path, di_freqIR)
    
        alpha_vis = a3 * ((ei_VIS[0] - 1.) / (ei_VIS[0] + 2.))
        alpha_scat = zeros((len(ei_SCAT)), dtype=complex)
        alpha_ir = zeros((len(ei_SCAT)), dtype=complex)
    
        for x in range(len(ei_SCAT)):
            alpha_scat[x] = a3 * ((ei_SCAT[x] - 1.) / (ei_SCAT[x] + 2.))
            alpha_ir[x]  = a3 * ((ei_IR[x] - 1.) / (ei_IR[x] + 2.))
        #E0 generated from SFG, E1 generated from visable incident, E2 generated from ir incident
        E0, E1, E2 = [zeros((len(ei_SCAT),3,3), dtype=complex) for x in range(3)]
        #E1 = zeros((3,3), dtype = complex)
        FG0, FG1, FG2 = [zeros((len(ei_SCAT),3,3,3),dtype=complex) for x in range(3)]
        #FG1 = zeros((3,3,3), dtype = complex)
        for n in range(len(ei_SCAT)):
            for l in range(3): # cycle over each pertubation direction
                e = zeros((3), dtype=complex)
                e[l] = 1.0
                # generate point dipole for nanoparticle
                mu0, mu1, mu2  = [zeros((3), dtype=complex) for x in range(3)]
                for i in range(3):
                    mu1[i] += alpha_vis * e[i]
                    mu0[i] += alpha_scat[n] * e[i]
                    mu2[i] += alpha_ir[n] *  e[i]
                # calculate E

                for i in range(3):
                    for j in range(3): # calculate e-field
                        E0[n][l][i] += 3.0 * r[i] * r[j] * mu0[j] / r5
                        if i==j: E0[n][l][i] -= mu0[j] / r3
                        E1[n][l][i] += 3.0 * r[i] * r[j] * mu1[j] / r5
                        if i==j: E1[n][l][i] -= mu1[j] / r3
                        E2[n][l][i] += 3.0 * r[i] * r[j] * mu2[j] / r5
                        if i==j: E2[n][l][i] -= mu2[j] / r3
                        #calculate FG
                        for k in range(3):
                            FG0[n][l][i][j] -= 15.0 * r[i]*r[j]*r[k]*mu0[k] / r7
                            FG1[n][l][i][j] -= 15.0 * r[i]*r[j]*r[k]*mu1[k] / r7
                            FG2[n][l][i][j] -= 15.0 * r[i]*r[j]*r[k]*mu2[k] / r7
                            if i==j:
                                FG0[n][l][i][j] += 3*r[k]*mu0[k]/r5
                                FG1[n][l][i][j] += 3*r[k]*mu1[k]/r5
                                FG2[n][l][i][j] += 3*r[k]*mu2[k]/r5
                            if i==k:
                                FG0[n][l][i][j] += 3*r[j]*mu0[k]/r5
                                FG1[n][l][i][j] += 3*r[j]*mu1[k]/r5
                                FG2[n][l][i][j] += 3*r[j]*mu2[k]/r5
                            if j==k:
                                FG0[n][l][i][j] += 3*r[i]*mu0[k]/r5
                                FG1[n][l][i][j] += 3*r[i]*mu1[k]/r5
                                FG2[n][l][i][j] += 3*r[i]*mu2[k]/r5

        return E0,E1,E2,FG0,FG1,FG2

    # generate alpha if none given
    #if alpha is None and SFG == False and hyperR == False:
    else:
        if alpha == None:
            if ei is None:
                freq = NM2HART(freq_VIS)
                ei = expdie('/usr/repo/applications/dim/dim/dielectrics/Ag_jc', freq)
#                ei = -1.985800+0.285400j # silver at 354nm
#                ei = -9.1297+0.3089j # silver at 488nm
            e0 = 1.0
            alpha = a3 * (( ei - e0 ) / ( ei + 2*e0 ))

        # calculate E and FG
        E   = 3. * einsum('i,j,jl->li', r, r, dt) * alpha / r5
        E  -= einsum('ij,jl->li', dt, dt) * alpha / r3
        FG  = -15. * einsum('i,j,k,kl->lij', r, r, r, dt) * alpha / r7
        FG += 3. * einsum('ij,k,kl->lij', dt, r, dt) * alpha / r5
        FG += 3. * einsum('ik,j,kl->lij', dt, r, dt) * alpha / r5
        FG += 3. * einsum('jk,i,kl->lij', dt, r, dt) * alpha / r5

        # return calculated E and FG
        return E, FG

def return_Raman_intensity(a, factor=1., **kwargs):
    '''Calculates the Raman intensity given a polarizability tensor derivative, a.'''
    from ..constants import BOHR2ANGSTROM as B2A
    from numpy import absolute, einsum
    cross1 = 0.; cross2 = 0.; diag = 0
    iso    = einsum('...aa->...', a) / 3.
    cross1 = einsum('...jk,...jk->...', a, a.conjugate())
    cross2 = einsum('...jk,...kj->...', a, a.conjugate())
    diag   = einsum('...jj,...kk->...', a, a.conjugate())
    ani = absolute((3./4.) * ( cross1 + cross2 ) - (1./2.) * diag)
    antisymm = 1.5 * einsum('...ik,...ik->...', a, a.conjugate())
    return (( 45. * iso * iso.conjugate() + 7. * ani + 5. * antisymm) / 45. ).real * B2A**4 * factor


def return_Raman_fixedframe(a, factor=1., incident='z', **kwargs):
    '''Calculates the Raman intensity given a fixed-frame regime where incident
       is the direction of the incident light. Assuming we are looking at back
       scattering for the intensity.'''
    from ..constants import BOHR2ANGSTROM as B2A
    from numpy import einsum
    if incident.lower() == 'z':
        S = ( a[0][0]*a[0][0].conjugate() + a[0][1]*a[0][1].conjugate()
            + a[1][0]*a[1][0].conjugate() + a[1][1]*a[1][1].conjugate() )
    elif incident.lower() == 'y':
        S = ( a[0][0]*a[0][0].conjugate() + a[2][0]*a[2][0].conjugate()
            + a[0][2]*a[0][2].conjugate() + a[2][2]*a[2][2].conjugate() )
    elif incident.lower() == 'x':
        S = ( a[1][1]*a[1][1].conjugate() + a[1][2]*a[1][2].conjugate()
            + a[2][1]*a[2][1].conjugate() + a[2][2]*a[2][2].conjugate() )
    return S * B2A**4 * factor


def calc_chi2_sfg(f, beta, **kwargs):
    '''
    Calculate the macroscopic chi2 tensor from the microscopic beta
    taken from J. Phys. Chem. B, 2004, 108 (11), 3548
    '''

    from numpy import where
    from numpy import zeros
    from math import radians, cos, sin

    # define chi2
    chi2 = zeros((f.nmodes,7),dtype='complex')

    rkwargs = return_kwargs
    theta = radians(rkwargs('theta',**kwargs))
    psi = radians(rkwargs('psi',**kwargs))

    for i in range(f.nmodes):

        # PPP - polarization
        #1 CHI2_ZZZ
        
        chi2[i,0] = (cos(theta)**3) * beta[i,2,2,2] \
                  + sin(theta) * sin(psi) * (beta[i,1,2,2] + beta[i,2,1,2] + beta[i,2,2,1]) \
                  - sin(theta) * cos(psi) * (beta[i,0,2,2] + beta[i,2,0,2] + beta[i,2,2,0])  \
                  + (sin(theta)**2) * cos(theta) * (sin(psi)**2) * (beta[i,1,1,2] + beta[i,1,2,1] + beta[i,2,1,1]) \
                  + (sin(theta)**2) * cos(theta) * (cos(psi)**2) * (beta[i,0,0,2] + beta[i,0,2,0] + beta[i,2,0,0])  \
                  - (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * (beta[i,0,1,2]+beta[i,0,2,1] + beta[i,1,0,2] + beta[i,1,2,0] + beta[i,2,0,1] + beta[i,2,1,0]) \
                  + (sin(theta)**3) * sin(psi) * (beta[i,0,0,1] + beta[i,0,1,0] + beta[i,1,0,0] - beta[i,1,2,2] - beta[i,2,1,2] - beta[i,2,2,1] ) \
                  + (sin(theta)**3) * cos(psi) * (beta[i,0,2,2] + beta[i,2,0,2] + beta[i,2,2,0] - beta[i,0,1,1] - beta[i,1,0,1] - beta[i,1,1,0] ) \
                  + (sin(theta)**3) * (sin(psi)**3) * (beta[i,1,1,1] - beta[i,0,0,1] - beta[i,0,1,0] - beta[i,1,0,0]) \
                  + (sin(theta)**3) * (cos(psi)**3) * (-beta[i,0,0,0] + beta[i,0,1,1] + beta[i,1,0,1] + beta[i,1,1,0]) 

        #2 CHI2_ZXX
        
        chi2[i,1] = ( (sin(theta)**2) * cos(theta)  * beta[i,2,2,2]
                  + cos(theta) * ( beta[i,2,0,0] + beta[i,2,1,1])
                  - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * (beta[i,1,1,2] + beta[i,1,2,1] + beta[i,2,1,1])
                  - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * (beta[i,0,0,2] + beta[i,0,2,0] + beta[i,2,0,0])
                  + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * (beta[i,0,1,2] + beta[i,0,2,1] + beta[i,1,0,2] + beta[i,1,2,0] + beta[i,2,0,1] + beta[i,2,1,0])
                  + sin(theta) * sin(psi) * (beta[i,1,1,1] + beta[i,1,0,0] - beta[i,2,1,2] - beta[i,2,2,1])
                  + sin(theta) * cos(psi) * (-beta[i,0,0,0] - beta[i,0,1,1] + beta[i,2,0,2] + beta[i,2,2,0])
                  + (sin(theta)**3) * sin(psi) * (-beta[i,0,0,1] - beta[i,0,1,0] - beta[i,1,0,0] + beta[i,1,2,2] + beta[i,2,1,2] + beta[i,2,2,1])
                  + (sin(theta)**3) * cos(psi) * ( beta[i,0,1,1] + beta[i,1,0,1] + beta[i,1,0,0] - beta[i,0,2,2] - beta[i,2,0,2] - beta[i,2,2,0])
                  + (sin(theta)**3) * (sin(psi)**3) * (-beta[i,1,1,1] + beta[i,0,0,1] + beta[i,0,1,0] + beta[i,1,0,0])
                  + (sin(theta)**3) * (cos(psi)**3) * (beta[i,0,0,0] - beta[i,0,1,1] - beta[i,1,0,1] - beta[i,1,1,0]) )

        chi2[i,1] = 0.5 * chi2[i,1]

        #3 CHI2_XZX
        
        chi2[i,2] = ( (sin(theta)**2) * cos(theta)  * beta[i,2,2,2]
                  + cos(theta) * ( beta[i,0,2,0] + beta[i,1,2,1])
                  - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * (beta[i,1,1,2] + beta[i,1,2,1] + beta[i,2,1,1])
                  - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * (beta[i,0,0,2] + beta[i,0,2,0] + beta[i,2,0,0])
                  + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * (beta[i,0,1,2] + beta[i,0,2,1] + beta[i,1,0,2] + beta[i,1,2,0] + beta[i,2,0,1] + beta[i,2,1,0])
                  + sin(theta) * sin(psi) * (beta[i,1,1,1] + beta[i,0,1,0] - beta[i,1,2,2] - beta[i,2,2,1])
                  + sin(theta) * cos(psi) * (-beta[i,0,0,0] - beta[i,1,0,1] + beta[i,0,2,2] + beta[i,2,2,0])
                  + (sin(theta)**3) * sin(psi) * (-beta[i,0,0,1] - beta[i,0,1,0] - beta[i,1,0,0] + beta[i,1,2,2] + beta[i,2,1,2] + beta[i,2,2,1])
                  + (sin(theta)**3) * cos(psi) * ( beta[i,0,1,1] + beta[i,1,0,1] + beta[i,1,1,0] - beta[i,0,2,2] - beta[i,2,0,2] - beta[i,2,2,0])
                  + (sin(theta)**3) * (sin(psi)**3) * (-beta[i,1,1,1] + beta[i,0,0,1] + beta[i,0,1,0] + beta[i,1,0,0])
                  + (sin(theta)**3) * (cos(psi)**3) * (beta[i,0,0,0] - beta[i,0,1,1] - beta[i,1,0,1] - beta[i,1,1,0]) )

        chi2[i,2] = 0.5 * chi2[i,2]


        #4 SSP - polarization
        # CHI2_XXZ SSP
        
        chi2[i,3] = ( (sin(theta)**2) * cos(theta)  * beta[i,2,2,2]
                  + cos(theta) * ( beta[i,0,0,2] + beta[i,1,1,2])
                  - (sin(theta)**2) * cos(theta) * (sin(psi)**2) * (beta[i,1,1,2] + beta[i,1,2,1] + beta[i,2,1,1])
                  - (sin(theta)**2) * cos(theta) * (cos(psi)**2) * (beta[i,0,0,2] + beta[i,0,2,0] + beta[i,2,0,0])
                  + (sin(theta)**2) * cos(theta) * sin(psi) * cos(psi) * (beta[i,0,1,2] + beta[i,0,2,1] + beta[i,1,0,2] + beta[i,1,2,0] + beta[i,2,0,1] + beta[i,2,1,0])
                  + sin(theta) * sin(psi) * (beta[i,1,1,1] + beta[i,0,0,1] - beta[i,1,2,2] - beta[i,2,1,2])
                  + sin(theta) * cos(psi) * (-beta[i,0,0,0] - beta[i,1,1,0] + beta[i,0,2,2] + beta[i,2,0,2])
                  + (sin(theta)**3) * sin(psi) * (-beta[i,0,0,1] - beta[i,0,1,0] - beta[i,1,0,0] + beta[i,1,2,2] + beta[i,2,1,2] + beta[i,2,2,1])
                  + (sin(theta)**3) * cos(psi) * ( beta[i,0,1,1] + beta[i,1,0,1] +beta[i,1,1,0] - beta[i,0,2,2] - beta[i,2,0,2] - beta[i,2,2,0])
                  + (sin(theta)**3) * (sin(psi)**3) * (-beta[i,1,1,1] + beta[i,0,0,1] + beta[i,0,1,0] + beta[i,1,0,0])
                  + (sin(theta)**3) * (cos(psi)**3) * (beta[i,0,0,0] - beta[i,0,1,1] - beta[i,1,0,1] - beta[i,1,1,0]) )

        chi2[i,3] = 0.5 * chi2[i,3]


        #5 CHI2_XYZ PSP
        
        chi2[i,4] = ( (cos(theta)**2) * (beta[i,0,1,2] - beta[i,1,0,2])
                  + (sin(theta)**2) * (sin(psi)**2) * ( beta[i,2,0,1] - beta[i,0,2,1])
                  + (sin(theta)**2) * (cos(psi)**2) * (beta[i,1,2,0] - beta[i,2,1,0] )
                  + (sin(theta)**2) * sin(psi) * cos(psi) * (beta[i,0,2,0] - beta[i,1,2,1] - beta[i,2,0,0] + beta[i,2,1,1])
                  + sin(theta) * cos(theta) * sin(psi) * (beta[i,0,1,1] - beta[i,0,2,2] - beta[i,1,0,1] + beta[i,2,0,2])
                  + sin(theta) * cos(theta) * cos(psi) * (-beta[i,0,1,0] + beta[i,1,0,0] - beta[i,1,2,2] + beta[i,2,1,2]) )

        chi2[i,4] = 0.5 * chi2[i,4]

        #6 CHI2_XZY PPS

        chi2[i,5] = ( (cos(theta)**2) * (beta[i,0,2,1] - beta[i,1,2,0])
                  + (sin(theta)**2) * (sin(psi)**2) * ( beta[i,2,1,0] - beta[i,0,1,2])
                  + (sin(theta)**2) * (cos(psi)**2) * (beta[i,1,0,2] - beta[i,2,0,1] )
                  + (sin(theta)**2) * sin(psi) * cos(psi) * (beta[i,0,0,2] - beta[i,1,1,2] - beta[i,2,0,0] + beta[i,2,1,1])
                  + sin(theta) * cos(theta) * sin(psi) * (beta[i,0,1,1] - beta[i,0,2,2] - beta[i,1,1,0] + beta[i,2,2,0])
                  + sin(theta) * cos(theta) * cos(psi) * (-beta[i,0,0,1] + beta[i,1,0,0] - beta[i,1,2,2] + beta[i,2,2,1]) )

        chi2[i,5] = 0.5 * chi2[i,5]
        

        #7 CHI2_ZXY PPS
       
        chi2[i,6] = ( (cos(theta)**2) * (beta[i,2,0,1] - beta[i,2,1,0])
                  + (sin(theta)**2) * (sin(psi)**2) * ( beta[i,1,2,0] - beta[i,1,0,2])
                  + (sin(theta)**2) * (cos(psi)**2) * (beta[i,0,1,2] - beta[i,0,2,1] )
                  + (sin(theta)**2) * sin(psi) * cos(psi) * (beta[i,0,0,2] - beta[i,0,2,0] - beta[i,1,1,2] + beta[i,1,2,1])
                  + sin(theta) * cos(theta) * sin(psi) * (beta[i,1,0,1] - beta[i,1,1,0] - beta[i,2,0,2] + beta[i,2,2,0])
                  + sin(theta) * cos(theta) * cos(psi) * (-beta[i,0,0,1] + beta[i,0,1,0] - beta[i,2,1,2] + beta[i,2,2,1]) )

        chi2[i,6] = 0.5 * chi2[i,6]
 
    
    chi1 = chi2
    chi2 = chi2 * chi2.conjugate() 
    return chi1, chi2


def print_datafile(f,chi2,**kwargs):
    '''Print the datafile to file or screen.'''
    from .functions import sum_lorentzian
    from numpy import where
    from numpy import linspace
    import matplotlib.pyplot as plt

    # Format strings
    #d  = '{0:>7.2f}'
    #t  = 'Tdip {0[0]:g} {0[1]:g} {0[2]:g}'
    #en = 'Energy {0:g}'

    # Set to render text with LaTeX, edit font properties, define figure
    #plt.rc('text', usetex=True)
    #plt.rc('font', **{'family':'serif', 'serif': ['Times'], 'size': 20})
    #fig = plt.figure()

    fwhm = 20
    # Define the points to plot, then plot
    domain = linspace(0, f.v_frequencies[-1]*1.5, num=2000)

    rkwargs = return_kwargs
    direction=rkwargs('direction',**kwargs)
    outfile=rkwargs('output',**kwargs)

    scale = 1
    # gamma == Full width at half max / 2
    gm = fwhm / 2
    y = sum_lorentzian(domain,f.v_frequencies, chi2[:,direction].real, gm)

    unNormy =y
#    graphoff = False
#    if not graphoff:
#        print (max(y))
#        #y = y/max(y)
#        sub = fig.add_subplot(111)
#        sub.plot(domain, y*scale, 'r')

        # Title, lables, and limits
        #sub.set_title('Title')
#        lab = r'SFG Intensity (a.u.)'
#        sub.set_ylabel(lab)
#        sub.set_ylim(0,1)
#        sub.set_xlabel(r'Wavenumber (cm$^{-1}$)')
        # plot to screen
        #plt.show()
    for i in range (len(domain)):
        print("{0}   {1} ".format(domain[i], unNormy[i]),file=outfile)
    outfile.close()

def expdie(name, freqs):
    '''Retrieves the experimental dielectric parameters for this type.
    If none were found, return complex zero's for each frequency'''
    from sys import argv
    from os.path import join, dirname, realpath
    from ..constants import HART2NM

    def loadtxt(filename):
        '''Loads reads a text file as a list of floats.
        Each column of data is a row in the returned data list.'''
        # Read in the entire file but skip comments. Split by columns
        with open(filename) as fl:
            f = [x.strip().split() for x in fl if x[0] != '#']
        # Now, reorganize from row x col ==> col x row
        return [[float(nums[n]) for nums in f] for n in xrange(len(f[0]))]

    die = []
    if name:
        try:
            expdata = loadtxt(name)
        except IOError:
            # If not found as given, look in the DIM library
            try:
                p = join(dirname(realpath(argv[0])), 'dielectrics', name)
                expdata = loadtxt(p)
            except IOError:
                exit('Cannot find dielectric file for '+name)

        # Spline the data
        realknots = spline(expdata[0], expdata[1])
        imagknots = spline(expdata[0], expdata[2])
        # For each frequency, interpolate
        try:
            for om in freqs:
                real = interpolate(expdata[0], expdata[1], realknots, HART2NM(om))
                imag = interpolate(expdata[0], expdata[2], imagknots, HART2NM(om))
                die.append(complex(real, imag))
        except TypeError:
            real = interpolate(expdata[0], expdata[1], realknots, HART2NM(freqs))
            imag = interpolate(expdata[0], expdata[2], imagknots, HART2NM(freqs))
            die = real + 1.j*imag
    else:
        die = [complex(0.0) for i in xrange(len(freqs))]

    return die


def spline(x, y):
    '''Splines the given data and returns the knots of the cubic spline.
    This spline is "natural", as in the bounary derivative is zero.'''


    #####################################################
    # First set up the three bands of the matrix to solve
    #####################################################

    # Get the length of the data
    n = len(x) - 1

    # Set up the tridiagonal equations
    c = [x[i] - x[i-1] for i in xrange(1,n)] + [0.0]
    c[0] = 0.0
    d = [1.0] + [2.0 * ( x[i+1] - x[i-1] ) for i in xrange(1,n)] + [1.0]
    e = [0.0] + [x[i+1] - x[i] for i in xrange(1,n)]

    # Set up the knots to solve for
    k = [0.0] + [6.0 * ( ( y[i+1] - y[i] ) / ( x[i+1] - x[i] ) 
                       - ( y[i] - y[i-1] ) / ( x[i] - x[i-1] ) )
                                                 for i in xrange(1,n)] + [0.0]

    ##########################################################
    # Solve the tridiagonal system to return the knots
    # This is like a pure python drop-in for dgtsv from LAPACK
    ##########################################################

    # Get the diagonal length
    n += 1

    # Decompose the matrix
    for i in xrange(1, n):
        lam    = c[i-1] / d[i-1]
        d[i]  -= lam * e[i-1]
        c[i-1] = lam

    # Use the decomposed matrix to solve for the knots
    for i in xrange(1,n):
        k[i] = k[i] - c[i-1] * k[i-1]
    k[n-1] = k[n-1] / d[n-1]
    for i in xrange(n-2, -1, -1):
        k[i] = ( k[i] - e[i] * k[i+1] ) / d[i]

    return k

def interpolate(x, y, knots, xvalue):
    '''Given the raw x and y data and the knots for that data, interpolate and
    return the y value associated with the requested x value'''

    ########################################################
    # Determine the index range between which our value lies
    ########################################################

    # The first two choices are for when the value matches our extremma
    if xvalue == x[0]:
        i = 0
    elif xvalue == x[-1]:
        i = len(x) - 2
    else:
        # Determine bounds
        ascending = x[-1] >= x[0]
        iLeft = -1
        iRight = len(x)

        # Find the index for the value immeadiately
        # below or equal to the requested value
        while iRight - iLeft > 1:
            # Compute a midpoint, and replace either the lower limit
            # or the upper limit, as appropriate.
            iMid = ( iLeft + iRight ) // 2
            if ascending == ( xvalue >= x[iMid] ):
                iLeft = iMid
            else:
                iRight = iMid

        # iLeft the index
        i = iLeft

    # Make sure the index falls in the correct window
    i = max(min(i, len(x)-2), 0) 
    
    ###################################
    # Interpolate using the given knots
    ###################################

    h = x[i+1] - x[i]
    a = ( x[i+1] - xvalue ) / h
    b = ( xvalue - x[i] )   / h
    return ( a * y[i] + b * y[i+1] 
         + ( ( a**3 - a ) * knots[i] + ( b**3 - b ) * knots[i+1] )
         * ( h**2 ) / 6.0 )

def symm(T):
    '''Retuns the symmetric and anti-symmetric parts
    of the second order tensor, T.'''
    from numpy import einsum
    Ts = 0.5 * ( T + einsum('...ab->...ba', T) )
    Ta = 0.5 * ( T - einsum('...ab->...ba', T) )
    return Ts, Ta

def check_kwargs(string, val=True, **kwargs):
    '''Check if a string is contained within the given kwargs
    and return True if its value is equal to val.'''
    if string in kwargs:
        return (kwargs[string] == val)
    else:
        return False

def return_kwargs(string, val=False, **kwargs):
    '''Check if a string is contained within the given kwargs
    and return the value of that kwargs, otherwise return val.'''
    if string in kwargs:
        return kwargs[string]
    else:
        return val

def plot_cars_raman(xc, yc, xr, yr, fwhm=20, minx=400, maxx=2000, step=1,
    fs=18, lw=1.5, **kwargs):

    from matplotlib import pyplot as plt
    from .functions import sum_lorentzian
    from numpy import arange, log10
    rkwargs = return_kwargs

    # Generate lorentzian for Raman
    x = arange(minx,maxx,step)
    yr = sum_lorentzian(x, xr, yr, fwhm=fwhm)

    # Scale spectra by appropriate factors
    scalec = int(log10(yc.max()))/3 * 3
    scaler = int(log10(yr.max()))/3 * 3
    yc = yc * 10**(-scalec)
    yr = yr * 10**(-scaler)

    # Plot CARS and Raman
    fig = plt.figure()

    sub1 = fig.add_subplot(111)
    sub1.axhline(linewidth=1, color='k')
    sub1.plot(x, yr,  'r-', lw=lw)
    sub1.set_ylabel('Raman ($10^{'+str(scaler)+'} cm^2/sr$)',
                    fontsize=fs, color='r')
    sub1.set_xlabel('wavenumber ($cm^{-1}$)', fontsize=fs)

    sub2 = sub1.twinx()
    sub2.plot(xc, yc,  'b-', lw=lw)
    sub2.set_ylabel('CARS ($10^{'+str(scalec)+'} cm^6 s^2 /{photon}^2 sr$)',
                    fontsize=fs, color='b')
    sub2.set_xlim(minx,maxx)

    # Are inverted axes needed?
    if rkwargs('invertx', True, **kwargs):
        plt.gca().invert_xaxis()
    if rkwargs('inverty', False, **kwargs):
        plt.gca().invert_yaxis()

    plt.show()
    return

def gen_field (f, radius=None, r=None, alpha=None, ei=None, e0=1., freq0=None, Einc=None, hyperpol=False, **kwargs):
    '''Generates the field and field gradient at the center of mass of
    the molecule (cm) given the the vector from the NP to the molecule (r)
    and the polarizability of the nanoparticle (alpha).
    e is the pertubating field vector, and is taken as 1 if no fields given.'''

    from ..constants import HART2NM, WAVENUM2HART, HART2EV
    from ..constants import KRONECKER3 as dt
    from numpy import array, zeros, einsum
    from math import sqrt

    gamma = return_kwargs('gamma', 50., **kwargs)
    if Einc is None: Einc = array([1.,1.,1.])

    # calculate r and a3
    if r is None: r = array([35.,0.,0.])
    #r1 = sqrt( ( r * r.conj() ).sum().real )
    # numpy-1.13 does not allow conjugate of a real number -- Pengchong
    r1 = sqrt( ( r * r ).sum() )
    r3 = r1 * r1 * r1
    r5 = r3 * r1 * r1
    r7 = r5 * r1 * r1
    r9 = r7 * r1 * r1
    if radius is None: radius = r1 # a is the radius of the sphere
    a3 = radius**3

    # generate alpha if none given
    # Zhongwei: take the 2\omega case into account for SEHRS
    # To make this consistent with the Camden group, for SEHRS, we pick the
    # scattered frequency to be on resonance with the Ag isotropic sphere
    #if freq0 is None: freq0 = f.e_frequencies[0]
    try: 
        if not f.b_e_frequencies[0]:
           if freq0 is None: freq0 = f.e_frequencies[0]
        else:
           if hyperpol: # incident frequency
              if freq0 is None: freq0 = f.b_e_frequencies[0]
           else: # scattered frequency 
              if freq0 is None: freq0 = f.e_frequencies[0]
    except IndexError:
        if freq0 is None: freq0 = f.e_frequencies[0]
    except TypeError:
        pass
    # Zhongwei: for the frequency dependent case, make sure the scattered
    # frequency is (about) twice larger than the incident frequency 
    # Pengchong: only do this for hyperpol; otherwise b_e_frequencies is an empty array, and 
    #   numpy-v1.13 does not allow index 0 for empty array
    if hyperpol:
        if (f.e_frequencies[0] < 2.1*f.b_e_frequencies[0]) and (f.e_frequencies[0] > 1.9*f.b_e_frequencies[0]):
           pass
        else:
           exit(
           '*****************************************************************************\n'
           'WARNING: For SEHRS simulations, the scattered frequency should be about twice\n'
           'larger than the incident frequency. Make sure to obtain the correct scattered\n'
           'frequency by changing the "group" value in the "collect_tensor_derivatives". \n'
           '*****************************************************************************')
    # Zhongwei: take the 2\omega case into account for SEHRS
    # To make this consistent with the Camden group, for SEHRS, we pick the
    # scattered frequency to be on resonance with the Ag isotropic sphere
    #if freq0 == 0.: freq0 = 0.1328 # If static, use 343 nm for Ag dielectric
    if freq0 == 0.:
       if hyperpol: # incident frequency
          freq0 = 0.0664 # If static and 2\omega, use 686 nm for Ag dielectric
          # alternatively
          #freq0 = 0.1328 # If static and 2\omega, use 343 nm for Ag dielectric
       else: #scattered frequency
          freq0 = 0.1328 # If static and \omega, use 343 nm for Ag dielectric
          # alternatively
          #freq0 = 0.2656 # If static and \omega, use 171 nm for Ag dielectric
    if alpha is None:
        if ei is None:
            if return_kwargs('gamma', **kwargs):
                ei = drude_dielectric(HART2EV(freq0), wp=HART2EV(freq0), gamma=gamma)
            else:
                #ei = expdie('/usr/repo/applications/dim/dim/dielectrics/Ag_jc', freq0)
                ei = expdie('/usr/repo/applications/dim/dim/dielectrics/Au_jc', freq0)
        alpha = [a3 * (( ei - e0 ) / ( ei + 2*e0 ))]

    # calculate E and FG
    # NB: traceless contributions not needed for a "flat" surface
    E   = 3. * einsum('i,l,a,l->ali', r, r, alpha, Einc) / r5
    if not return_kwargs('flatsurface', **kwargs):
        E  -= einsum('il,a,l->ali', dt, alpha, Einc) / r3
    FG  = -15. * einsum('i,j,l,a,l->alij', r, r, r, alpha, Einc) / r7
    if not return_kwargs('flatsurface', **kwargs):
        FG += 3. * einsum('ij,l,a,l->alij', dt, r, alpha, Einc) / r5
        FG += 3. * einsum('il,j,a,l->alij', dt, r, alpha, Einc) / r5
        FG += 3. * einsum('jl,i,a,l->alij', dt, r, alpha, Einc) / r5

    # return calculated E and FG
    return E[0], FG[0]

def drude_dielectric(w, wp=3.54, gamma=0.5):
    '''Calculates the dielectric assuming a drude function.
    w, wp and gamma are given in eV.'''
    ei = 2. - 4. * ( wp**2 / ( w**2 + 1j*w*gamma ))
    return ei

def lorentzian_normal_sum(x,y,z,dx,dy,dz,mr,mi,bxr,byr,bzr,bxi,byi,bzi):
    from numpy import zeros,sqrt,pi
    fxr=(0.5*bxr)/((x-dx)**2+(0.5*bxr)**2)/pi
    fxi=(0.5*bxi)/((x-dx)**2+(0.5*bxi)**2)/pi
    fyr=(0.5*byr)/((y-dy)**2+(0.5*byr)**2)/pi
    fyi=(0.5*byi)/((y-dy)**2+(0.5*byi)**2)/pi
    fzr=(0.5*bzr)/((z-dz)**2+(0.5*bzr)**2)/pi
    fzi=(0.5*bzi)/((z-dz)**2+(0.5*bzi)**2)/pi
    tx,ty,tz,t = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    tx.real,tx.imag = fxr,fxi
    ty.real,ty.imag = fyr,fyi
    tz.real,tz.imag = fzr,fzi
    t.real = tx.real*ty.real*tz.real*mr
    t.imag = tx.imag*ty.imag*tz.imag*mi
    field=zeros((len(x),3,3),dtype=complex)
    field[:,2,2]=t
    fdxr = -bxr*(x-dx)/((x-dx)**2+(0.5*bxr)**2)**2/pi
    fdxi = -bxi*(x-dx)/((x-dx)**2+(0.5*bxi)**2)**2/pi
    fdyr = -byr*(y-dy)/((y-dy)**2+(0.5*byr)**2)**2/pi
    fdyi = -byi*(y-dy)/((y-dy)**2+(0.5*byi)**2)**2/pi
    fdzr = -bzr*(z-dz)/((z-dz)**2+(0.5*bzr)**2)**2/pi
    fdzi = -bzi*(z-dz)/((z-dz)**2+(0.5*bzi)**2)**2/pi
    tdx,tdy,tdz = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    ldx,ldy,ldz = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    tdx.real,tdx.imag = fdxr,fdxi
    tdy.real,tdy.imag = fdyr,fdyi
    tdz.real,tdz.imag = fdzr,fdzi
    ldx.real = tdx.real*ty.real*tz.real*mr
    ldx.imag = tdx.imag*ty.imag*tz.imag*mi
    ldy.real = tdy.real*tx.real*tz.real*mr
    ldy.imag = tdy.imag*tx.imag*tz.imag*mi
    ldz.real = tdz.real*tx.real*ty.real*mr
    ldz.imag = tdz.imag*tx.imag*ty.imag*mi
    fg=zeros((len(x),3,3,3),dtype=complex)
    fg[:,2,2,0]=ldx
    fg[:,2,2,1]=ldy
    fg[:,2,2,2]=ldz
    return field, fg

def lorentzian_normal_mag(x,y,z,dx,dy,dz,mr,mi,bxr,byr,bzr,bxi,byi,bzi):
    from numpy import zeros,sqrt
    fxr=(0.5*bxr)**2/((x-dx)**2+(0.5*bxr)**2)
    fxi=(0.5*bxi)**2/((x-dx)**2+(0.5*bxi)**2)
    fyr=(0.5*byr)**2/((y-dy)**2+(0.5*byr)**2)
    fyi=(0.5*byi)**2/((y-dy)**2+(0.5*byi)**2)
    fzr=(0.5*bzr)**2/((z-dz)**2+(0.5*bzr)**2)
    fzi=(0.5*bzi)**2/((z-dz)**2+(0.5*bzi)**2)
    tx,ty,tz,t = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    tx.real,tx.imag = fxr,fxi
    ty.real,ty.imag = fyr,fyi
    tz.real,tz.imag = fzr,fzi
    t.real = tx.real*ty.real*tz.real*mr
    t.imag = tx.imag*ty.imag*tz.imag*mi
    field=zeros((len(x),3,3),dtype=complex)
    field[:,2,2]=t
    fdxr = -bxr**2*(x-dx)*0.5/((x-dx)**2+(0.5*bxr)**2)**2
    fdxi = -bxi**2*(x-dx)*0.5/((x-dx)**2+(0.5*bxi)**2)**2
    fdyr = -byr**2*(y-dy)*0.5/((y-dy)**2+(0.5*byr)**2)**2
    fdyi = -byi**2*(y-dy)*0.5/((y-dy)**2+(0.5*byi)**2)**2
    fdzr = -bzr**2*(z-dz)*0.5/((z-dz)**2+(0.5*bzr)**2)**2
    fdzi = -bzi**2*(z-dz)*0.5/((z-dz)**2+(0.5*bzi)**2)**2
    tdx,tdy,tdz = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    ldx,ldy,ldz = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    tdx.real,tdx.imag = fdxr,fdxi
    tdy.real,tdy.imag = fdyr,fdyi
    tdz.real,tdz.imag = fdzr,fdzi
    ldx.real = tdx.real*ty.real*tz.real*mr
    ldx.imag = tdx.imag*ty.imag*tz.imag*mi
    ldy.real = tdy.real*tx.real*tz.real*mr
    ldy.imag = tdy.imag*tx.imag*tz.imag*mi
    ldz.real = tdz.real*tx.real*ty.real*mr
    ldz.imag = tdz.imag*tx.imag*ty.imag*mi
    fg=zeros((len(x),3,3,3),dtype=complex)
    #fg[:,2,2,0],fg[:,2,0,2]=ldx,ldx
    #fg[:,2,2,1],fg[:,2,2,1]=ldy,ldy
    fg[:,2,2,0]=ldx
    fg[:,2,2,1]=ldy
    fg[:,2,2,2]=ldz
    return field, fg

def gaussian(x,y,z,dx,dy,dz,mr,bxr,byr,bzr,bxi,byi,bzi):
    from numpy import zeros,exp,sqrt
    frx=exp(-(x-dx)**2*0.5/bxr**2)
    fix=exp(-(x-dx)**2*0.5/bxi**2)
    fry=exp(-(y-dy)**2*0.5/byr**2)
    fiy=exp(-(y-dy)**2*0.5/byi**2)
    frz=exp(-(z-dz)**2*0.5/bzr**2)
    fiz=exp(-(z-dz)**2*0.5/bzi**2)
    tx,ty,tz,t = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    tx.real,tx.imag = frx,fix
    ty.real,ty.imag = fry,fiy
    tz.real,tz.imag = frz,fiz
    t.real = tx.real*ty.real*tz.real*mr
    t.imag = tx.imag*ty.imag*tz.imag*mi
    field=zeros((len(x),3,3),dtype=complex)
    field[:,2,2]=t
    #first derivative
    frdx = -(x-dx)/bxr**2*frx
    fidx = -(x-dx)/bxi**2*fix
    frdy = -(y-dy)/byr**2*fry
    fidy = -(y-dy)/byi**2*fiy
    frdz = -(z-dz)/bzr**2*frz
    fidz = -(z-dz)/bzi**2*fiz
    ldx,ldy,ldz = zeros(len(x),dtype=complex),zeros(len(x),dtype=complex),zeros(len(x),dtype=complex)
    ldx.real = tdx.real*ty.real*tz.real*mr
    ldx.imag = tdx.imag*ty.imag*tz.imag*mi
    ldy.real = tdy.real*tx.real*tz.real*mr
    ldy.imag = tdy.imag*tx.imag*tz.imag*mi
    ldz.real = tdz.real*tx.real*ty.real*mr
    ldz.imag = tdz.imag*tx.imag*ty.imag*mi
    fg=zeros((len(x),3,3,3),dtype=complex)
    fg[:,2,2,0]=ldx
    fg[:,2,2,1]=ldy
    fg[:,2,2,2]=ldz
    return field, fg
