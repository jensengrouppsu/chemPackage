#! /usr/bin/env python

from __future__ import print_function

def dressed_spectroscopy(obj, r=None, debug=False, gradient=False, magnetic=False, spectrum='Raman',
    ws=False, lplot=True,  **kwargs):
    '''
    Calculates the specified spectrum from a given ChemData object "obj".
    NB. ChemData object obj must have the polarzability derivaties already
    calculated. These can be obtained through numerical differentiation with
    "nmodes2numdiff.py" followed by the function obj.collect_roa_derivatives().
    Spectrum type may be specified by the *spectrum* key, defaulted to Raman.

    Addtional options are:

    E         ==> Explicitly gives the local electric field (a 3x3 tensor).
    FG        ==> Explicitly gives the local electric field-gradient (a 3x3x3 tensor).
    r         ==> Generates the fields from an isotropic sphere, where the vector
                  r describes the separation of the molecule from the center of the
                  sphere and |r| is the radius of the sphere.
    gradient  ==> Whether to use field-gradient effects (default False).
    magnetic  ==> Whether to use magnetic field effects (default False)
    debug     ==> Turn on print statements.
    spectrum  ==> Specify type of spectrum to be plotted (may be one of
                  "Raman" - default, "ROA", "CID", "CARS", or "hyperRaman").
    noplot    ==> Returns the frequency and SERS, SEROA or SECARS values as arrays.
    sticks    ==> Plots the sticks spectrum along with the lorenztian broadened
                  spectrum (defaul True). If specified along with *noplot* it returns
                  the intensity per mode, if set to False it returns the lorenztian
                  spectrum.
    Atensor   ==> Whether to use the dressed-Atensor contributions to the ROA/CID
                  spectrum (defaul True).
    Gtensor   ==> Whether to use the dressed-Gtensor contributions to the ROA/CID
                  spectrum (default True).
    unrestricted ==> Use the unrestricted Raman, ROA, or CID expressions. This makes
                  a distinction between the Roman-type and script-type A- and G-tensors.
    ws        ==> Whether to not use the approximation omega_L .neq. omega_S.
                  May only be used for fields from an isotropic sphere (and
                  therefore, r MUST be specified). For ROA or CID, may only be used
                  with the *unrestricted* key.
    lplot     ==> logical variable that determines the format of the output. If true
                   this function will return a plot of the requested spectrum, if false
                   this function will return the x,y values for the spectrum instead. 
                   Suggested to always use false because plotting is sacred and is best 
                   done yourself, just sayin' :/

    The required tensor derivatives are:
    Spectrum    Normal / E-field only        Field-gradient
    ----------  ---------------------        --------------
    Raman       alpha                        alpha, A, C
    ROA         alpha, A, G                  alpha, A, C, G, D
    CID         alpha, A, G                  alpha, A, C, G, D
    CARS        alpha                        alpha, A, C
    hyperRaman  beta                         beta, ed_ed_eq, ed_eq_ed,
                                             eq_ed_ed, ed_eq_eq, eq_ed_eq,
                                             eq_eq_ed, eq_eq_eq
    '''

    from numpy import array, zeros, einsum
    from ..constants import ANGSTROM2BOHR as A2B, KRONECKER3 as dt
    from .dressed_func import plot, generate_field, print_field, raman_cross_section
    from .dressed_func import check_kwargs as ckwargs, return_kwargs as rkwargs
    from .dressed_func import gen_field, hyperraman_cross_section
    from copy import deepcopy
    from .translate_tensors import translate_polarizabilities

    # Define some initial valutes
    if ws and not rkwargs('unrestricted', **kwargs):
        kwargs.update({'unrestricted': True})
    spectrum = spectrum.upper()
    try:
        R = rkwargs('com', array(obj.center_of_nuc_charge*A2B,dtype=float), **kwargs)
    except AttributeError:
        R = zeros((3))

    # Unfold A, C and D tensors into 3x3x3 or 3x3x3x3 tensors
    obj.fold_Atensor('UNFOLD')
    obj.fold_Ctensor('UNFOLD')
    obj.fold_Dtensor('UNFOLD')

    # Translate polarizabilities if needed
#    if (abs(R)>1e-2).any(): obj = translate_polarizabilities(obj, R)

    # Generate Fields
    E = rkwargs('E', zeros((3,3)), **kwargs)
    B = rkwargs('B', zeros((3,3)), **kwargs)
    FG = rkwargs('FG', zeros((3,3,3)), **kwargs)
    Es = rkwargs('Es', None, **kwargs)
    Bs = rkwargs('Bs', None, **kwargs)
    FGs = rkwargs('FGs', None, **kwargs)
    # Zhongwei: field and field-gradient for 2\omega
    E2 = rkwargs('E2', zeros((3,3)), **kwargs)
    FG2 = rkwargs('FG2', zeros((3,3,3)), **kwargs)
    if r is not None:
        # Treat the hyperraman case to be consistent with the Camden group
        if (spectrum=='HYPERRAMAN' or spectrum=='HSERS' or spectrum=='HRS'):
           E, FG = gen_field(obj, r=r, freqs=[obj.v_frequencies if ws else None][0], hyperpol=True, **kwargs)
           E2, FG2 = gen_field(obj, r=r, freqs=[obj.v_frequencies if ws else None][0], **kwargs)
           #if ws: Es, FGs = gen_field(obj, r=r, **kwargs)
        else:
           E, FG = gen_field(obj, r=r, freqs=[obj.v_frequencies if ws else None][0], **kwargs)
           if ws: Es, FGs = gen_field(obj, r=r, **kwargs)
    # Print values used for local electric field and gradient
    if debug and not ((E==0.).all() and ((FG==0.).all() or not gradient)) and not ws:
        if gradient:
            print_field(E, FG)
        else:
            print_field(E)
    # Add delta to generate enhancement factor
    E = dt + E
    if ws and Es is not None: Es = dt + Es
    #Zhongwei: also do this for E2
    E2 = dt + E2

    # Build dressed-tensors
    if (spectrum=='RAMAN' or spectrum=='NRS' or spectrum=='SERS'
        or spectrum=='CARS' or spectrum=='SECARS'):
        print( "hello my dudes")
        if magnetic:
            aD = dress_alpha_plus_mag(obj, E, B, FG, gradient=gradient, magnetic=magnetic, Es=Es, Bs=Bs, FGs=FGs)
        else:
            aD = dress_alpha(obj, E, FG, gradient=gradient, Es=Es, FGs=FGs)
    if (spectrum=='HYPERRAMAN' or spectrum=='HSERS' or spectrum=='HRS'):
        #Zhongwei: include E2 and FG2 for hyper-SERS 
        BD = dress_hyperRaman_beta(obj, E, FG, E2, FG2, gradient=gradient)
    elif not (spectrum=='RAMAN' or spectrum=='NRS' or spectrum=='SERS'
            or spectrum=='CARS' or spectrum=='SECARS'):
        AD = dress_A(obj, E, FG, gradient=gradient)
        GD = dress_G(obj, E, FG, gradient=gradient)
    else:
        AD = zeros((len(aD),3,3,3), dtype=aD.dtype)
        GD = zeros((len(aD),3,3), dtype=aD.dtype)

    # Dress script-type tensors if needed (only works with unrestricted formalism)
    if rkwargs('unrestricted', **kwargs) and (spectrum=='ROA' or spectrum=='SEROA'
       or spectrum=='CID' or spectrum=='SECID'):
        if r is not None and ws:
            AsD = dress_A(obj, Es, FGs, gradient=gradient)
            GsD = -einsum('iab->iba', dress_G(obj, Es, FGs, gradient=gradient))
        else:
            AsD = dress_A(obj, E, FG, gradient=gradient)
            GsD = -einsum('iab->iba', dress_G(obj, E, FG, gradient=gradient))
    else:
        AsD = None
        GsD = None

    # If we don't need certain tensors, set to zero
    if ckwargs('Atensor', False, **kwargs): AD = zeros((len(aD),3,3,3))
    if ckwargs('Gtensor', False, **kwargs): GD = zeros((len(aD),3,3))

    # Return intensity
    if (spectrum=='HYPERRAMAN' or spectrum=='HSERS' or spectrum=='HRS'):
        # Zhongwei: calculate either static or frequency-dependent hyper-SERS intensities
        intensity = hyperraman_cross_section(obj.v_frequencies, BD, b_e_frequencies=obj.b_e_frequencies, **kwargs)
    elif spectrum=='CARS' or spectrum=='SECARS':
        return CARS_intensity(obj, aD, **kwargs)
    else:
        if rkwargs('unrestricted', **kwargs):
            intensity, sers = unrestricted_ROA(aD, AD, GD, obj.e_frequencies, AsD, GsD,
                                               debug=debug, spectrum=spectrum, **kwargs)
        else:
            print("maybe i'm here?")
            intensity, sers = ROA_intensity(aD, AD, GD, obj.e_frequencies, debug=debug,
                                            spectrum=spectrum, **kwargs)

        # Convert intensity to cross-section if not CIDs
        sers = raman_cross_section(obj.v_frequencies, sers, e_frequencies=obj.e_frequencies, **kwargs)
        if spectrum!='CID' and spectrum!='SECID':
            intensity = raman_cross_section(obj.v_frequencies, sers, e_frequencies=obj.e_frequencies, **kwargs)

    # Print intensitiy to standard output
    if debug:
        maxi = max(abs(intensity.max()), abs(intensity.min()))
        maxs = abs(sers.max())
        temps = 'ROA'
        if ckwargs('spectrum', 'CID', **kwargs): temps = 'CID'
        print ('')
        print (' Frequency (cm-1)     (SE)'+temps+'           (SE/N)RS')
        for i in range(obj.nmodes):
            if (abs(intensity[i])>(maxi/10.)): #or (abs(sers[i])>(maxs/10.)):
                print (' {0:10.2f}          {1: 10.6e}      {2:10.6e}'.format(
                       obj.v_frequencies[i], intensity[i], sers[i]))

    # Choosing the format of the output.
    if not lplot:
        # return the x,y values of the spectrum
        if spectrum=='RAMAN' or spectrum=='NRS' or spectrum=='SERS':
            return obj.v_frequencies, sers
        else:
            return obj.v_frequencies, intensity
    else:
        # return the plot of the spectum. i.e. you're doing it wrong. (no offense)
        if spectrum=='RAMAN' or spectrum=='NRS' or spectrum=='SERS':
            return plot(obj.v_frequencies, sers, spectrum=spectrum, **kwargs)
        else:
            return plot(obj.v_frequencies, intensity, spectrum=spectrum, **kwargs)

def ROA_intensity(a, A, G, omega, debug=False, geom='BACKWARD', component='all', **kwargs):
    '''ROA intensties as defined per J. Chem. Phys., 127, 134101 (2007)
    and page 347 of the Barron (2004) text.'''

    from ..constants import LEVICIVITA3 as eps, LIGHT_AU as c, BOHR2ANGSTROM as B2A
    import sys
    from .dressed_func import check_kwargs as ckwargs
    from numpy import einsum

    # Zhongwei: in case we only want certain components for Raman/SERS
    if component != 'all':
       if len(component) == 2:
          print('Considering only ' + component.upper() + ' component')
          sym_change = { 'X' : '0', 'Y' : '1', 'Z' : '2', }
          sym_index = []
          for i in range(len(component)):
              sym_index.append(int(sym_change[component[i].upper()]))
          for i in range(len(a)):
              for j in range(3):
                  for k in range(3):
                      if j == sym_index[0] and k == sym_index[1]:
                         pass
                      else:
                         a[i,j,k] = 0.0
       else:
          print('The component ' + component.upper() + ' is invalid for polarizability.')
          print('Please use the XY format.')
          import sys
          sys.exit()

    geom = geom.upper()
    k = (1./45.) * B2A**4

    # Set some strings
    ndim = a.ndim
    smode = ''
    if ndim == 3: smode = 'i'
    string1 = '{0}aa,{0}bb->{0}'.format(smode)
    string2 = '{0}ab,{0}ab->{0}'.format(smode)
    string3 = '{0},{0}ab,acd,{0}cdb->{0}'.format(smode)

    # Generate invariants
    a2  = (1./9.) * ( einsum(string1, a, a.conjugate()) ).real
    Ba  = 1.5 * ( einsum(string2, a, a.conjugate())  ).real
    Ba -= 0.5 * ( einsum(string1, a, a.conjugate())  ).real
    BG  = 1.5 * ( 1.j * einsum(string2, a, G.conjugate()) ).imag
    BG -= 0.5 * ( 1.j * einsum(string1, a, G.conjugate()) ).imag
    aG  = (1./9.) * ( 1.j * einsum(string1, a, G.conjugate()) ).imag
    BA  = 0.5 * ( einsum(string3, omega, a, eps, A.conjugate()) ).real

    # Return ROA intensity based on given scattering geometry
    if geom == 'FORWARD':
        roa   = (4*k/c) * ( 45*aG + BG - BA )
        raman = k * ( 45*a2 + 7*Ba )
    elif geom == 'NINETY_X':
        roa   = (k/c) * ( 45*aG + 7*BG + BA )
        raman = k * ( 45*a2 + 7*Ba )
    elif geom == 'NINETY_Z':
        roa   = (6*k/c) * ( BG - BA/3 )
        raman = 6*k*Ba
    elif geom == 'BACKWARD':
        roa   = (24*k/c) * ( BG + BA/3 )
        raman = k * ( 45*a2 + 7*Ba )
    else:
        print ('Unrecognized geometry specification.')
        print ('Must be one of:')
        print ('FORWARD, NINETY_X, NINETY_Z, BACKWARD')
        sys.exit()

    # Return intensities
    if ckwargs('spectrum', 'CID', **kwargs):
        return roa/raman, raman
    else:
        return roa, raman

def unrestricted_ROA(a, A, G, omega, A_s=None, G_s=None, debug=False,
    angle=180, scattering='ICP', polarization='U', **kwargs):
    '''ROA intensities as defined in Nafie, Vibrational Optical Activity, 2011. pp 133-167.'''

    from copy import deepcopy
    from .dressed_func import symm, check_kwargs as ckwargs
    from numpy import einsum
    from ..constants import LEVICIVITA3 as eps, LIGHT_AU as c
    from ..constants import BOHR2ANGSTROM as B2A
    import sys

    # Assume script-type tensors are identical to Roman-type if not present
    if A_s is None: A_s = deepcopy(A)
    if G_s is None: G_s = -einsum('...ab->...ba', G)

    # Define factor k
    k = (1./180.) * B2A**4

    # Generate a temp string
    st = ''
    st2 = ''
    if a.ndim==3 and A.ndim==4 and G.ndim==3:
        st = 'i'
        st2 = '->i'

    # Get the 2nd rank tensor contraction of epsilon and A
    A16 = einsum('acd,{0}cdb->{0}ab'.format(st), eps, A)
    A12 = einsum('abc,{0}dcd->{0}ab'.format(st), eps, A)
    A16_s = einsum('acd,{0}cdb->{0}ab'.format(st), eps, A_s)
    A12_s = einsum('abc,{0}dcd->{0}ab'.format(st), eps, A_s)

    # Return symmetric and anti-symmetric tensors
    aas, aaa = symm(a)
    Gs, Ga = symm(G)
    G_ss, G_sa = symm(G_s)
    eA1s, eA1a = symm(A16)
    eA1_ss, eA1_sa = symm(A16_s)
    eA2s, eA2a = symm(A12)
    eA2_ss, eA2_sa = symm(A12_s)

    # Generate tensor invariants
    a2 = (einsum('{0}aa,{0}bb{1}'.format(st,st2), aas, aas.conjugate())).real / 9.
    Bsa = 1.5 * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aas, aas.conjugate())).real
    Bsa -= 0.5 * (einsum('{0}aa,{0}bb{1}'.format(st,st2), aaa, aaa.conjugate())).real
    Baa = 1.5 * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aaa, aaa.conjugate())).real
    aG = ( 1j * einsum('{0}aa,{0}bb{1}'.format(st,st2), aas, Gs.conjugate())).imag / 9.
    BsG = 1.5 * ( 1j * einsum('{0}ab,{0}ab{1}'.format(st,st2), aas, Gs.conjugate())).imag
    BsG -= 0.5 * ( 1j * einsum('{0}aa,{0}bb{1}'.format(st,st2), aas, Gs.conjugate())).imag
    BaG = 1.5 * ( 1j * einsum('{0}ab,{0}ab{1}'.format(st,st2), aaa, Ga.conjugate())).imag
    BsA = 0.5 * omega * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aas, eA1s.conjugate())).real
    BaA = 0.5 * omega * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aaa, eA1a.conjugate())).real
    BaA += 0.5 * omega * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aaa, eA2a.conjugate())).real
    aG_s = ( 1j * einsum('{0}aa,{0}bb{1}'.format(st,st2), aas, G_ss.conjugate())).imag / 9.
    BsG_s = 1.5 * ( 1j * einsum('{0}ab,{0}ab{1}'.format(st,st2), aas, G_ss.conjugate())).imag
    BsG_s -= 0.5 * ( 1j * einsum('{0}aa,{0}bb{1}'.format(st,st2), aas, G_ss.conjugate())).imag
    BaG_s = 1.5 * ( 1j * einsum('{0}ab,{0}ab{1}'.format(st,st2), aaa, G_sa.conjugate())).imag
    BsA_s = 0.5 * omega * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aas, eA1_ss.conjugate())).real
    BaA_s = 0.5 * omega * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aaa, eA1_sa.conjugate())).real
    BaA_s += 0.5 * omega * (einsum('{0}ab,{0}ab{1}'.format(st,st2), aaa, eA2_sa.conjugate())).real

    # Generate intensities based on geometry

    # FORWARD
    if (angle==0 and scattering=='ICP' and polarization=='U'):
        raman = 2*k * ( 90*a2 + 14*Bsa + 10*Baa )
        roa = (4*k/c) * ( 90*aG + 14*BsG + 10*BaG + 2*BsA - 2*BaA - 90*aG_s + 10*BsG_s - 10*BaG_s - 6*BsA_s + 2*BaA_s )

    elif (angle==0 and scattering=='SCP' and polarization=='U'):
        raman = 2*k * ( 90*a2 + 14*Bsa + 10*Baa )
        roa = (4*k/c) * ( 90*aG - 10*BsG + 10*BaG - 6*BsA - 2*BaA - 90*aG_s - 14*BsG_s - 10*BaG_s + 2*BsA_s + 2*BaA_s )

    elif (angle==0 and scattering=='DCP1'):
        raman = 2*k * ( 90*a2 +  2*Bsa + 10*Baa )
        roa = (4*k/c) * ( 90*aG +  2*BsG + 10*BaG - 2*BsA - 2*BaA - 90*aG_s -  2*BsG_s - 10*BaG_s - 2*BsA_s + 2*BaA_s )

    elif (angle==0 and scattering=='DCP2'):
        raman = 2*k * (  0*a2 + 12*Bsa +  0*Baa )
        roa = (4*k/c) * (  0*aG + 12*BsG +  0*BaG + 4*BsA + 0*BaA +  0*aG_s + 12*BsG_s +  0*BaG_s - 4*BsA_s + 0*BaA_s )

    # RIGHT ANGLE
    elif (angle==90 and scattering=='ICP' and polarization=='P'):
        raman = 2*k * ( 45*a2 + 7*Bsa + 5*Baa )
        roa = (4*k/c) * ( 45*aG + 7*BsG + 5*BaG + 1*BsA - 1*BaA + 0*aG_s + 0*BsG_s + 0*BaG_s + 0*BsA_s + 0*BaA_s )

    elif (angle==90 and scattering=='ICP' and polarization=='D'):
        raman = 2*k * ( 0*a2 + 6*Bsa + 10*Baa )
        roa = (4*k/c) * ( 0*aG + 6*BsG + 10*BaG - 2*BsA + 2*BaA + 0*aG_s + 0*BsG_s + 0*BaG_s + 0*BsA_s + 0*BaA_s )

    elif (angle==90 and scattering=='DCP1'):
        raman = k * ( 45*a2 + 13*Bsa + 15*Baa )
        roa = (2*k/c) * ( 45*aG + 13*BsG + 15*BaG - 1*BsA + 1*BaA - 45*aG_s - 13*BsG_s - 15*BaG_s - 1*BsA_s - 1*BaA_s )

    elif (angle==90 and scattering=='DCP2'):
        raman = k * ( 45*a2 + 13*Bsa + 15*Baa )
        roa = (2*k/c) * ( 45*aG + 13*BsG + 15*BaG - 1*BsA + 1*BaA + 45*aG_s + 13*BsG_s + 15*BaG_s + 1*BsA_s + 1*BaA_s )

    # BACKWARD
    elif (angle==180 and scattering=='SCP' and polarization=='U'):
        raman = 2*k * ( 90*a2 + 14*Bsa + 10*Baa )
        roa = (4*k/c) * ( - 90*aG + 10*BsG - 10*BaG + 6*BsA + 2*BaA - 90*aG_s - 14*BsG_s - 10*BaG_s + 2*BsA_s + 2*BaA_s )

    elif (angle==180 and scattering=='DCP1'):
        raman = 2*k * (  0*a2 + 12*Bsa +  0*Baa )
        roa = (4*k/c) * (  0*aG + 12*BsG +  0*BaG + 4*BsA + 0*BaA +  0*aG_s - 12*BsG_s +  0*BaG_s + 4*BsA_s + 0*BaA_s )

    elif (angle==180 and scattering=='DCP2'):
        raman = 2*k * ( 90*a2 +  2*Bsa + 10*Baa )
        roa = (4*k/c) * ( 90*aG +  2*BsG + 10*BaG - 2*BsA - 2*BaA + 90*aG_s +  2*BsG_s + 10*BaG_s + 2*BsA_s - 2*BaA_s )

    elif (angle==180 and scattering=='ICP' and polarization=='U'):
        raman = 2*k * ( 90*a2 + 14*Bsa + 10*Baa )
        roa = (4*k/c) * ( 90*aG + 14*BsG + 10*BaG + 2*BsA - 2*BaA + 90*aG_s - 10*BsG_s + 10*BaG_s + 6*BsA_s - 2*BaA_s )

    else:
        print ('Unrecognized option for "angle", "scattering" and "polarization".')
        print ('Exiting.')
        sys.exit()

    # Return intensities
    if ckwargs('spectrum', 'CID', **kwargs):
        return roa/raman, raman
    else:
        return roa, raman

def CARS_intensity(obj, aD, lifetime=0.00005, noplot=False, **kwargs):
    '''
    Generate the CARS tensor and intensity.
    '''

    from numpy import array, zeros, arange, einsum
    from ..constants import WAVENUM2HART, NM2HART
    from .dressed_func import cars_cross_section
    from .dressed_func import return_Raman_intensity, raman_cross_section
    from .dressed_func import return_kwargs as rkwargs, plot_cars_raman

    # Automatically readjust lifetime to match fwhm if given
    lifetime = rkwargs('fwhm', 20, **kwargs) / 400000.

    # Get number of input frequencies and initialize arrays
    minx = rkwargs('minx', 400., **kwargs)
    maxx = rkwargs('maxx', 2000., **kwargs)
    step = rkwargs('step', 1., **kwargs)
    om2 = arange(minx, maxx+step, step)
    om2 = NM2HART(obj.e_frequencies[0] - WAVENUM2HART(om2))

    # Calculate some frequencies
    wko = WAVENUM2HART(obj.v_frequencies)
    w1  = obj.e_frequencies[0]
    w2  = NM2HART(om2)
    w12 = w1 - w2
    freq = array(w12 / WAVENUM2HART, dtype=float)

    # Generate Y
    Y1  = einsum('iab,icd->iabcd', aD, aD)
    Y2  = einsum('iad,icb->iabcd', aD, aD)
    d1  = 1. / ( einsum('w,m->wm', -w12, zeros((len(wko)))+1) + wko + 1j*lifetime )
    d2  = 1. / ( einsum('w,m->wm',  w12, zeros((len(wko)))+1) + wko + 1j*lifetime )
    Y   = 2. * ( einsum('iabcd,wi->wiabcd', Y1, d1 )
               + einsum('iabcd,wi->wiabcd', Y2, d2 ))

    # Calculate isotropic Y
    Yiso = (1./15.) * ( einsum('wiabba->wi', Y) + einsum('wiaa->wi', einsum('wiabcc->wiab', Y))
         + einsum('wiabab->wi', Y) )
    S2 = ( Yiso * Yiso.conjugate() ).real
    S2 = einsum('wi->w', S2)

    # Convert to cross-section
    S2 = cars_cross_section(freq, S2, obj.e_frequencies)

    # Return intensities if not plotting
    if noplot: return freq, S2

    # Collect Raman intensity and cross-section
    raman = return_Raman_intensity(aD)
    raman = raman_cross_section(obj.v_frequencies, raman, obj.e_frequencies)

    # Plot
    return plot_cars_raman(freq, S2, obj.v_frequencies, raman, **kwargs)

def dress_alpha(obj, E, FG, gradient=False, Es=None, FGs=None):
    '''Dress the polarzability tensor with E and FG.'''
    from numpy import einsum, zeros
    from sys import exit
    from copy import deepcopy
    if Es is None: Es = deepcopy(E)
    if FGs is None: FGs = deepcopy(FG)
    st1 = ''
    if E.ndim == 3 and FG.ndim == 4: st1 = 'i'
    st2 = ''
    if Es.ndim==3 and FGs.ndim==4: st2 = 'i'
    aD  = einsum('{1}ac,icd,{0}bd->iab'.format(st1,st2), Es, obj.polarizability, E)
    #########################################################################
    # Zhongwei: print out the \alpha components values for analyzing the
    # orientation information from SERS
    #for i in range(len(aD)):
    #    print ('Mode:', obj.v_frequencies[i])
    #    for j in range(3):
    #        for k in range(3):
    #            print (j,k,aD[i][j][k].real,aD[i][j][k].imag)
    #    print ('-----------------------------------------------------------')
    #########################################################################
    if gradient:
        try:
            aD += (1./3.) * einsum('{1}ac,icde,{0}bde->iab'.format(st1,st2), Es, obj.atensor, FG)
            aD += (1./3.) * einsum('{1}acd,iecd,{0}be->iab'.format(st1,st2), FGs, obj.atensor, E)
            aD += (1./9.) * einsum('{1}acd,icdeg,{0}beg->iab'.format(st1,st2), FGs, obj.ctensor, FG)
        except ValueError:
            exit('ERROR: A- and C-tensors are needed for field-gradient effects!')
    return aD

def dress_alpha_plus_mag(obj, E, B, FG, gradient=False, magnetic=False, Es=None, Bs=None, FGs=None):
    '''Dress the polarzability tensor with E and FG.'''
    from numpy import einsum, zeros
    from sys import exit
    from copy import deepcopy
    if Es is None: Es = deepcopy(E)
    if Bs is None: Bs = deepcopy(B)
    if FGs is None: FGs = deepcopy(FG)
    st1 = ''
    if E.ndim == 3 and FG.ndim == 4: st1 = 'i'
    st2 = ''
    if Es.ndim==3 and FGs.ndim==4: st2 = 'i'
    print( st1, st2, E.ndim, FG.ndim, obj.polarizability.ndim)
    #print( obj.polarizability)
    aD  = einsum('{1}ac,icd,{0}bd->iab'.format(st1,st2), Es, obj.polarizability, E)
    #########################################################################
    # Zhongwei: print out the \alpha components values for analyzing the
    # orientation information from SERS
    #for i in range(len(aD)):
    #    print ('Mode:', obj.v_frequencies[i])
    #    for j in range(3):
    #        for k in range(3):
    #            print (j,k,aD[i][j][k].real,aD[i][j][k].imag)
    #    print ('-----------------------------------------------------------')
    #########################################################################
    if gradient:
        try:
            aD += (1./3.) * einsum('{1}ac,icde,{0}bde->iab'.format(st1,st2), Es, obj.atensor, FG)
            aD += (1./3.) * einsum('{1}acd,iecd,{0}be->iab'.format(st1,st2), FGs, obj.atensor, E)
            aD += (1./9.) * einsum('{1}acd,icdeg,{0}beg->iab'.format(st1,st2), FGs, obj.ctensor, FG)
        except ValueError:
            exit('ERROR: A- and C-tensors are needed for field-gradient effects!')
    if magnetic:
        try:
            aD +=  einsum('{1}ac,icd,{0}bd->iab'.format(st1,st2), Es, obj.gtensor, B)
            aD +=  einsum('{1}ac,icd,{0}bd->iab'.format(st1,st2), Bs, obj.gtensor, E)
        except ValueError:
            exit('ERROR: G-tensor is needed for magnetic field effects!')

    return aD



def dress_A(obj, E, FG, gradient=False):
    '''Dress the A-tensor with E and FG.'''
    from numpy import einsum
    from sys import exit
    st = ''
    if E.ndim==3 and FG.ndim==4: st = 'i'
    try:
        AD  = einsum('{0}ad,idbc->iabc'.format(st), E, obj.atensor)
    except ValueError:
        exit ('ERROR: A-tensor not found!')
    if gradient:
        try:
            AD += (1./3.) * einsum('{0}ade,idebc->iabc'.format(st), FG, obj.ctensor)
        except ValueError:
            exit ('ERROR: C-tensor needed for field-gradient effects!')
    return AD

def dress_G(obj, E, FG, gradient=False):
    '''Dress the G-tensor with E and FG.'''
    from numpy import einsum
    from sys import exit
    st = ''
    if E.ndim==3 and FG.ndim==4: st = 'i'
    try:
        GD  = einsum('{0}ac,icb->iab'.format(st), E, obj.gtensor)
    except ValueError:
        exit ('ERROR: G-tensor not found!')
    if gradient:
        try:
            GD += (1./3.) * einsum('{0}acd,icdb->iab'.format(st), FG, obj.dtensor)
        except ValueError:
            exit ('ERROR: D-tensor needed for field-gradient effects!')
    return GD

def dress_hyperRaman_beta(obj, E, FG, E2, FG2, gradient=False):
    '''Dress the beta tensor with E and FG for hyperSERS.'''
    from numpy import einsum, array
    from sys import exit
    st = ''
    if E.ndim==3 and FG.ndim==4: st = 'i'
    st2 = ''
    if E2.ndim==3 and FG2.ndim==4: st2 = 'i'
    try:
        # Zhongwei: make the order consistent with the Raman case
        #           1. scattering frequency   (\omega_s --> E2)
        #           2. response property      (\beta --> hyperpolarizaiblity)
        #           3. 1st incident frequency (\omega_1 --> E)
        #           4. 2nd incident frequency (\omega_2 --> E)
        #           Note: we treat \omega_s = \omega_1 + \omega_2 = 2\omega
        BD = einsum('{0}ad,idef,{0}be,{1}cf->iabc'.format(st2,st), E2, obj.hyperpolarizability, E, E)
    except ValueError:
        exit ('ERROR: beta tensor not found!')

    if gradient:
        # Expand quadrupole indices
        iQ = [[0, 1, 2], [1, 3, 4], [2, 4, 5]]
        DDQ = array([[[[[obj.beta_ddq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in range(3)]
                    for j in range(3)] for i in range(obj.nmodes)])
        DQD = array([[[[[obj.beta_dqd[i][j][k][l] for l in range(3)] for k in iQ[m]] for m in range(3)]
                    for j in range(3)] for i in range(obj.nmodes)])
        QDD = array([[[[[obj.beta_qdd[i][j][k][l] for l in range(3)] for k in range(3)] for j in iQ[m]]
                    for m in range(3)] for i in range(obj.nmodes)])
        DQQ = array([[[[[[obj.beta_dqq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in iQ[n]]
                    for n in range(3)] for j in range(3)] for i in range(obj.nmodes)])
        QDQ = array([[[[[[obj.beta_qdq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in range(3)]
                    for j in iQ[n]] for n in range(3)] for i in range(obj.nmodes)])
        QQD = array([[[[[[obj.beta_qqd[i][j][k][l] for l in range(3)] for k in iQ[m]] for m in range(3)]
                    for j in iQ[n]] for n in range(3)] for i in range(obj.nmodes)])
        QQQ = array([[[[[[[obj.beta_qqq[i][j][k][l] for l in iQ[m]] for m in range(3)] for k in iQ[n]]
                    for n in range(3)] for j in iQ[o]] for o in range(3)] for i in range(obj.nmodes)])

        BD += (1./3.) * einsum('{0}ad,idefg,{0}be,{1}cfg->iabc'.format(st2,st), E2, DDQ, E, FG)
        BD += (1./3.) * einsum('{0}ad,idefg,{0}bef,{1}cg->iabc'.format(st2,st), E2, DQD, FG, E)
        BD += (1./3.) * einsum('{0}ade,idefg,{0}bf,{1}cg->iabc'.format(st2,st), FG2, QDD, E, E)
        BD += (1./9.) * einsum('{0}ad,idefgh,{0}bef,{1}cgh->iabc'.format(st2,st), E2, DQQ, FG, FG)
        BD += (1./9.) * einsum('{0}ade,idefgh,{0}bf,{1}cgh->iabc'.format(st2,st), FG2, QDQ, E, FG)
        BD += (1./9.) * einsum('{0}ade,idefgh,{0}bfg,{1}ch->iabc'.format(st2,st), FG2, QQD, FG, E)
        BD += (1./27.) * einsum('{0}ade,idefghj,{0}bfg,{1}chj->iabc'.format(st2,st), FG2, QQQ, FG, FG)
    return BD
