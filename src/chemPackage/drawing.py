from __future__ import print_function
#from drawinginfo import __doc__

def drawAtoms(chemObj):
    ''' Draw all the atoms in 3D using VPython

    .. note::
        Requires that VPython and the Microsoft Default Font Package
        is also installed.

        '''
    import vis as vi
    from .constants import atomic_color
    # Loop over all the atoms in the system, plotting a sphere for each one
    # Radius and color are determined by dictionary lookups
    for i in xrange(chemObj.nallatoms):
        ele = chemObj.allatoms[i]
        rad = chemObj.allradii('vis')[i]
        vi.sphere(pos=chemObj.allcoords[i], radius=rad,
                    color=atomic_color(ele))
    # Draw axes to give the user a frame of reference
    vi.curve(pos=[(0,0,-15),(0,0,15)], radius=.2, color=vi.color.white)
    vi.text(pos=(0,0,-15), text='-Z', height=2, color=vi.color.white)
    vi.text(pos=(0,0,15), text='+Z', height=2, color=vi.color.white)
    vi.curve(pos=[(-15,0,0),(15,0,0)], radius=.1, color=vi.color.blue)
    vi.text(pos=(-15,0,0), text='-X', height=2, color=vi.color.blue)
    vi.text(pos=(15,0,0), text='+X', height=2, color=vi.color.blue)
    vi.curve(pos=[(0,-15,0),(0,15,0)], radius=.1, color=vi.color.green)
    vi.text(pos=(0,-15,0), text='-Y', height=2, color=vi.color.green)
    vi.text(pos=(0,15,0), text='+Y', height=2, color=vi.color.green)

def cubeFile (chemObj, xpara, ypara, zpara, file, angstrom=True,
              screen=False, smear=1.0, dir=4, sq=False, mag=False,
              type='static scattered', freq=0, logscale=False, qm=False,vec=False):
    ''' Generate a Gaussian cube file for use in 3D visualization programs

    *xpara*, *ypara*, and *zpara* define the 3D box you wish to draw.  Each of
    these should be 2 parameters: a min and a max.

    *vsize* is the size of the volume space element used in the grid, also
    known as a voxel.  The larger the voxel, the fewer the grid points resulting
    in a faster calculation but with lower resolution.  Users are currently
    limited to no more than 10,000,000 voxels to prevent generating cube files
    larger than 500MB.

    *angstrom* indicates that the grid values given are in Angstroms.  If
    set to False, we are assuming they have been given in atomic units (Bohr).

    The calculation type must be given and depending on the *type*, additional
    arguments must also be passed.  The following is a list of type allowed:

    ==================  ============================================
    ------------------  --------------------------------------------
    'ground image'      | Requires no other arguments
    'ES image'          | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)
                        | Requires frequency argument *freq*
                          which is the index of the desired frequency
                          in :py:attr:`~.e_frequencies`
                          Defaults to 0 (lowest frequency)

    'static scattered'  | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)

    'FD scattered'      | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)
                        | Requires frequency argument *freq*
                          which is the index of the desired frequency
                          in :py:attr:`~.e_frequencies`
                          Defaults to 0 (lowest frequency)
    ==================  ============================================


    *vec* will cause this method to return the two coordinate grids and a
    grid of e-field vectors instead of plotting them internally.  This is
    useful for manipulating the data outside the scope of this method.

    *min* sets the minimum to the scale bar.  Defaults to 0.

    *max* sets the maximum to the scale bar.  Defaults to the maximum of the
    e-field magnitudes^4

    *step* sets the setp size for the scale bar.  Defaults to 1.

    *sq* will use the electric field squared if set to True instead of E^4.
    Defaults to False.  This is mutually exclusive with *mag*.

    *mag* will use the magnitude of the e-field (no squaring, no 4th power) if
    set to True.  Defaults to False.  This is mutually exclusive with *sq*.

    *screen* will screen the distance between points when calculating the
    electric field using a gaussian smear.  Default is False (no screening).

    *smear* is the distance (in Bohr) that the grid points are smeared by when
    using screening.  Defaults to 1.0.
    '''
    from scipy.special import erf
    from .constants import BOHR2ANGSTROM, ANGSTROM2BOHR, atomic_radius, atomic_color
    from sys import exit
    from numpy import array, transpose, arange, empty_like, zeros,zeros_like, sqrt, log10, linspace
    #
    # Convert parameters to float arrays
    xpara = array(xpara, dtype=float)
    ypara = array(ypara, dtype=float)
    zpara = array(zpara, dtype=float)
    #
    # Set up screening value
    scrn = 0
    if screen:
        if 'DDA' in chemObj.key:
            scrn = 2
        else:
            scrn = 1
    if mag and sq:
        exit("Can only specify E^4, E^2, or E")
    #
    if not qm:
        # Convert to Bohr and transpose coordinates
        coords = transpose(chemObj.dim_coordinates*ANGSTROM2BOHR)
        print("Type: ", type)
        if dir < 4:
            print("Direction: ", ['X', 'Y', 'Z'][dir])

        #
        # Collect the coordinates and dipoles depending on the type of DIM
        # calculation chemObj is.
        if type == 'ground image':
            dipole = chemObj.dim_dipoles[type]
            if chemObj.dim_charges is not None: charges = chemObj.dim_charges[type]
        else:
            if dir == 4:
                exit("Need direction of polarization")
            if type == 'static scattered':
                dipole = chemObj.dim_dipoles[type][dir]
                if chemObj.dim_charges is not None: charges = chemObj.dim_charges[type][dir]
            else:
                print("Frequency: ", chemObj.e_frequencies[0])
                dipole = chemObj.dim_dipoles[type][freq][dir]
                if chemObj.dim_charges is not None: charges = chemObj.dim_charges[type][freq][dir]

        # if there are no DIM charges (ie a PIM job), the charges will be set to zero
        if chemObj.dim_charges is None: charges = zeros(len(dipole[:,0]))

        # Transpose the dipole array so we can do matrix-vector multiplication
        dipole = transpose(dipole)
        charges = transpose(charges)
        atom_radii = chemObj.dim_radii()*ANGSTROM2BOHR
        natoms=chemObj.ndim
        atmNum = chemObj.dim_atomic_numbers
    else:
    # With qm==True, calculate fields from hirshfeld partitioned polarizabilities/charges
    # -- Pengchong Oct. 2016
        coords = transpose(chemObj.coordinates*ANGSTROM2BOHR)
        atom_radii=chemObj.radii()*ANGSTROM2BOHR
        natoms = chemObj.natoms
        atmNum = chemObj.atomic_numbers
        dip=chemObj.hirshfeld_induced_dipoles_loc
        dipole=transpose(dip[:,dir,:])    
        charges = transpose(chemObj.hirshfeld_induced_charges[:,2])
    # Determine grid dimensions
    nx = round((xpara[1]-xpara[0])/xpara[2] + 1.)
    ny = round((ypara[1]-ypara[0])/ypara[2] + 1.)
    nz = round((zpara[1]-zpara[0])/zpara[2] + 1.)

#    if angstrom:
#        xpara = xpara*ANGSTROM2BOHR
#        ypara = ypara*ANGSTROM2BOHR
#        zpara = zpara*ANGSTROM2BOHR
    linex = linspace(xpara[0], xpara[1], nx)
    liney = linspace(ypara[0], ypara[1], ny)
    linez = linspace(zpara[0], zpara[1], nz)

    if angstrom:
        linex = ANGSTROM2BOHR(linex)
        liney = ANGSTROM2BOHR(liney)
        linez = ANGSTROM2BOHR(linez)
    numX = len(linex)
    numY = len(liney)
    numZ = len(linez)
    xc = ((xpara[1] - xpara[0])/2.0 + xpara[0])
    yc = ((ypara[1] - ypara[0])/2.0 + ypara[0])
    zc = ((zpara[1] - zpara[0])/2.0 + zpara[0])
    zc = ((zpara[1] - zpara[0])/2.0 + zpara[0])
    tV = numX * numY * numZ
#    print("Number of voxels: {0}".format(tV))
    if tV > 10000000:
        print("Too many voxels.  Increase their size to reduce their number.")
        exit(1)

##xing
    if vec:
        # Create a 1D array for each component of the electric field
        # This a the grid flattened
        x      = empty_like(linez)
        y      = empty_like(linez)
        efieldx = empty_like(linez)
        efieldx = array(efieldx, dtype='complex')
        efieldy = empty_like(linez)
        efieldy = array(efieldy, dtype='complex')
        efieldz = empty_like(linez)
        efieldz = array(efieldz, dtype='complex')

    else: 
        from f2py import drawing_math as dm
        x      = empty_like(linez)
        y      = empty_like(linez)
        efield = empty_like(linez)
        efield = array(efield, dtype='complex', order='f')

    if vec:
        from f2py import drawing_math_vectors as dmv
        from numpy import array,append
        efx,efy,efz=array([]),array([]),array([])
        for ix in xrange(numX):   
            x[:] = linex[ix]
            for iy in xrange(numY):
                y[:] = liney[iy]
                dmv(x, y, linez, coords, dipole, atom_radii, efieldx, efieldy, efieldz, charges, False, dir, scrn, smear)
                for iz in xrange(numZ):
                    efx=append(efx,efieldx[iz])
                    efy=append(efy,efieldy[iz])
                    efz=append(efz,efieldz[iz])
        from numpy import column_stack
        efield = column_stack((efx,efy,efz))
        return efield

    else:
        #
        # Print header to the cube file
        f = open(file, 'w')
        print('DRAWING CUBE FILE', file=f)
        print('OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z', file=f)
        print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(natoms, xpara[0], ypara[0], zpara[0]), file=f)
        print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(numX, xpara[2], 0.0, 0.0), file=f)
        print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(numY, 0.0, ypara[2], 0.0), file=f)
        print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(numZ, 0.0, 0.0, zpara[2]), file=f)
    #    atmNum = chemObj.dim_atomic_numbers
        for i in xrange(natoms):
            print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}'.format(atmNum[i], 0.0, coords[0][i], coords[1][i], coords[2][i]), file=f)
        #
        # Generate voxel information
        t = 0
        ef=[]
        for ix in xrange(numX):
            x[:] = linex[ix]
            for iy in xrange(numY):
                y[:] = liney[iy]
                dm(x, y, linez, coords, dipole, atom_radii, efield, charges, False, dir, scrn, smear)
                if mag:
                    efield = sqrt(efield)
                elif not sq:
                    efield = efield*efield
                if (logscale): efield = log10(efield)
                for iz in xrange(numZ):
                    print('{0:< 13.5E}'.format(efield[iz].real), end='', file=f)
                    ef.append(efield[iz])
                    if (t % 6 == 5):
                        print('\n', end='', file=f)
                    t += 1
        print('\n', end='', file=f)
        f.close()
        return ef

def drawField(chemObj, xpara, ypara, zpara, scale=1, calctype='static scattered',
                dir=4, max=None, min=0.0, step=1, smear=1.0, vec=False, freq=0,
                sq=False, mag=False, screen=False, draw=True, qm=False, cg=False):
    ''' Calculates the electric field on a 2D plane in space and displays
    a contour plot of the result

    .. note::
        Requires that matplotlib is installed

    One of the parameters must be a single value defining the plane to be viewed.
    The other two parameters are length 3 defined as:

    .. hlist::
        :columns: 3

        - min
        - max
        - stepsize

    These are used to define the dimensions of the plane.  Distances entered
    should be in angstroms.

    *scale* will scale the radius of the DIM atoms to make them their effective
    radius.  The default value is calculated according to
    J. Phys. Chem. C, 2008, 112 (40), 15697; Eq. 15  (0.7978845608028654)

    The calculation type must be given and depending on the *calctype*, additional
    arguments must also be passed.  The following is a list of type allowed:

    ==================  ============================================
    ------------------  --------------------------------------------
    'ground image'      | Requires no other arguments
    'ES image'          | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)
                        | Requires frequency argument *freq*
                          which is the index of the desired frequency
                          in :py:attr:`~.e_frequencies`
                          Defaults to 0 (lowest frequency)

    'static scattered'  | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)

    'FD scattered'      | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)
                        | Requires frequency argument *freq*
                          which is the index of the desired frequency
                          in :py:attr:`~.e_frequencies`
                          Defaults to 0 (lowest frequency)
    ==================  ============================================


    *vec* will cause this method to return the two coordinate grids and a
    grid of e-field vectors instead of plotting them internally.  This is
    useful for manipulating the data outside the scope of this method.

    *min* sets the minimum to the scale bar.  Defaults to 0.

    *max* sets the maximum to the scale bar.  Defaults to the maximum of the
    e-field magnitudes^4

    *step* sets the setp size for the scale bar.  Defaults to 1.

    *sq* will use the electric field squared if set to True instead of E^4.
    Defaults to False.  This is mutually exclusive with *mag*.

    *mag* will use the magnitude of the e-field (no squaring, no 4th power) if
    set to True.  Defaults to False.  This is mutually exclusive with *sq*.

    *screen* will screen the distance between points when calculating the
    electric field using a gaussian smear.  Default is False (no screening).

    *smear* is the distance (in Bohr) that the grid points are smeared by when
    using screening.  Defaults to 1.0.

    *draw* will draw the electric fields using matplotlib.  If set to False,
    it will instead return the grid points and efields at each grid point.

    '''
    from numpy import array, transpose, arange, zeros, zeros_like, reshape, log10, sqrt, shape, real
    from matplotlib import ticker
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from .constants import BOHR2ANGSTROM, ANGSTROM2BOHR, atomic_radius, atomic_color
    from sys import exit
    #
    # Convert parameters to float arrays
    xpara = array(xpara, dtype='float32')
    ypara = array(ypara, dtype='float32')
    zpara = array(zpara, dtype='float32')
    
    # Set up screening value
    scrn = 0
    if screen:
        if 'DDA' in chemObj.key:
            scrn = 2
        else:
            scrn = 1

    if not qm:    
        # Convert to Bohr and transpose coordinates
        coords = transpose(chemObj.dim_coordinates*ANGSTROM2BOHR)
        if cg:
            atom_radii = scale*chemObj.dim_cgradii*ANGSTROM2BOHR
        else:
            atom_radii = scale*chemObj.dim_radii()*ANGSTROM2BOHR

        print("Type: ", calctype)
        if dir < 4:
            try:
                print("Direction: ", ['X', 'Y', 'Z'][dir])
            except:
                exit("Must give valid direction index (0, 1, or 2)")

        # for retardation dir should match pol_vec of output file, if not, exit.
        print( chemObj.key['A_VEC'])
        if "RETARDATION" in chemObj.key:
            if dir != chemObj.key["POL_VEC"]:
                exit("Direction supplied must match polarizaton direction specified in calculation, the default pol_vec direction is y")

        if mag and sq:
            exit("Can only specify E^4, E^2, or E")

        #
        # Collect the coordinates and dipoles depending on the type of DIM
        # calculation chemObj is.
        if calctype == 'ground image':
            dipole = chemObj.dim_dipoles[calctype]
            if chemObj.dim_charges is not None: charges = chemObj.dim_charges[calctype]
        else:
            if dir == 4:
                exit("Need direction of polarization")
            if calctype == 'static scattered':
                dipole = chemObj.dim_dipoles[calctype][dir]
                if chemObj.dim_charges is not None: charges = chemObj.dim_charges[calctype][dir]
            else:
                print("Frequency: ", chemObj.e_frequencies[freq])
                dipole = chemObj.dim_dipoles[calctype][freq][dir]
                if chemObj.dim_charges is not None: charges = chemObj.dim_charges[calctype][freq][dir]

        # if there are no DIM charges (ie a PIM job), the charges will be set to zero
        if chemObj.dim_charges is None: charges = zeros(len(dipole[:,0]))

        # Transpose the dipole array so we can do matrix-vector multiplication
        dipole = transpose(dipole)
        charges = transpose(charges)
        ### QM
    else:
    # With qm==True, calculate fields from hirshfeld partitioned polarizabilities/charges
    # -- Pengchong Oct. 2016
        coords = transpose(chemObj.coordinates*ANGSTROM2BOHR)
        atom_radii=chemObj.radii()
        dip_tot=(chemObj.hirshfeld_induced_dipoles_loc+chemObj.hirshfeld_induced_dipoles_loc)
        dipole=transpose(dip_tot[:,dir,:])    
        charges = zeros(len(dipole[0,:]))

    # Determine which plane we are drawing by checking which parameter is only
    # length 1
    if xpara.size == 1:
        y, z, x, A, B = __build2DGrid(ypara, zpara, xpara)
        title = "X = " + str(xpara)
        labelA = "Y Axis"
        labelB = "Z Axis"
    elif ypara.size == 1:
        x, z, y, A, B = __build2DGrid(xpara, zpara, ypara)
        title = "Y = " + str(ypara)
        labelA = "X Axis"
        labelB = "Z Axis"
    elif zpara.size == 1:
        x, y, z, A, B = __build2DGrid(xpara, ypara, zpara)
        title = "Z = " + str(zpara)
        labelA = "X Axis"
        labelB = "Y Axis"
    else:
        exit("Must define a plane; only parameters can have a range")
    # Initialize the results array to zero and cast is as complex
    if vec:
        # Create a 1D array for each component of the electric field
        # This a the grid flattened
        efieldx = zeros_like(x)
        efieldy = zeros_like(x)
        efieldz = zeros_like(x)
        #tefield = zeros_like(x)
        efieldx = array(efieldx, dtype='complex')
        efieldy = array(efieldy, dtype='complex')
        efieldz = array(efieldz, dtype='complex')
        #tefield = array(efieldz, dtype='complex',order='f')
    else:
        efield = zeros_like(x)
        efield = array(efield, dtype='complex', order='f')
    nsolv = chemObj.key["NSOLV"]
    print('nsolv: ', nsolv)
    if vec:
        from .f2py import drawing_math_vectors as dmv
        from .f2py import drawing_math_vectors_ret as dmv_ret
        # Returns efieldx, efieldy, efieldz
        if "RETARDATION" in chemObj.key:
            avec = chemObj.key['A_VEC']
            dmv_ret(x, y, z, coords, dipole, atom_radii, efieldx, efieldy, efieldz,
                       charges, False, dir, scrn, smear, chemObj.e_frequencies[freq], nsolv, avec)
        else:
            dmv(x, y, z, coords, dipole, atom_radii, efieldx, efieldy, efieldz,
                 charges, False, dir, scrn, smear)
        # Stack and reshape into grid of vectors
        from numpy import column_stack
        efield = column_stack((efieldx, efieldy, efieldz))
        efield = reshape(efield, (A.shape[0], A.shape[1],3))
        # Return coordinate arrays and the vector array
        return A, B, efield

    from .f2py import drawing_math as dm
    from .f2py import drawing_math_ret as dm_ret
    
    if "RETARDATION" in chemObj.key:
        avec = chemObj.key['A_VEC']
        dm_ret(x, y, z,  coords,  dipole, atom_radii,efield , charges,
                False, dir, scrn, smear,chemObj.e_frequencies[freq], nsolv, avec) 
    else:
        print(len(x), len(y), len(z), coords.shape, dipole.shape, len(atom_radii),len(efield), len(charges)) 
        dm(x, y, z, coords, dipole, atom_radii, efield, charges, False, dir, scrn, smear)
    efield = reshape(efield, A.shape)

    if mag:
        efield = sqrt(efield)
    elif not sq:
        efield = efield*efield
    #
    # Creates the contour plot
    if draw:
        plt.figure(facecolor='#7F7F7F')
        #plt.axes(axisbg='k')
        pallete = cm.jet
        if max is None:
            cont = plt.contourf(A, B, log10(efield),
                                locator=ticker.MaxNLocator(500),
                                cmap=pallete,
                                )
        else:
            cont = plt.contourf(A, B, log10(efield), levels=arange(min,max,step, dtype=float),
                                extend='both',
                                cmap=pallete,
                                )
        CB = plt.colorbar(cont)
        if sq:
            CB.set_label('Electric Field Magnitude^2 (log(Atomic Units))')
        else:
            CB.set_label('Electric Field Magnitude^4 (log(Atomic Units))')
        plt.title(title + " (Angstroms)")
        plt.xlabel(labelA + " (Angstroms)")
        plt.ylabel(labelB + " (Angstroms)")
        plt.show()
        plt.close()
        del A, B, efield
        import gc
        gc.collect()
    else:
        return A, B, efield

def drawMagField(chemObj, xpara, ypara, zpara, scale=1, calctype='static scattered',
                dir=4, max=None, min=0.0, step=1, smear=1.0, vec=False, freq=0,
                sq=False, mag=False, screen=False, draw=True, qm=False, cg=False):
    ''' Calculates the magnetic field on a 2D plane in space and displays
    a contour plot of the result

    .. note::
        Requires that matplotlib is installed

    One of the parameters must be a single value defining the plane to be viewed.
    The other two parameters are length 3 defined as:

    .. hlist::
        :columns: 3

        - min
        - max
        - stepsize

    These are used to define the dimensions of the plane.  Distances entered
    should be in angstroms.

    *scale* will scale the radius of the DIM atoms to make them their effective
    radius.  The default value is calculated according to
    J. Phys. Chem. C, 2008, 112 (40), 15697; Eq. 15  (0.7978845608028654)

    The calculation type must be given and depending on the *calctype*, additional
    arguments must also be passed.  The following is a list of type allowed:

    ==================  ============================================
    ------------------  --------------------------------------------
    'ground image'      | Requires no other arguments
    'ES image'          | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)
                        | Requires frequency argument *freq*
                          which is the index of the desired frequency
                          in :py:attr:`~.e_frequencies`
                          Defaults to 0 (lowest frequency)

    'static scattered'  | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)

    'FD scattered'      | Requires direction argument *dir*
                          where x=0, y=1, and z=2
                          Defaults to 0 (X direction)
                        | Requires frequency argument *freq*
                          which is the index of the desired frequency
                          in :py:attr:`~.e_frequencies`
                          Defaults to 0 (lowest frequency)
    ==================  ============================================


    *vec* will cause this method to return the two coordinate grids and a
    grid of e-field vectors instead of plotting them internally.  This is
    useful for manipulating the data outside the scope of this method.

    *min* sets the minimum to the scale bar.  Defaults to 0.

    *max* sets the maximum to the scale bar.  Defaults to the maximum of the
    e-field magnitudes^4

    *step* sets the setp size for the scale bar.  Defaults to 1.

    *sq* will use the electric field squared if set to True instead of E^4.
    Defaults to False.  This is mutually exclusive with *mag*.

    *mag* will use the magnitude of the e-field (no squaring, no 4th power) if
    set to True.  Defaults to False.  This is mutually exclusive with *sq*.

    *screen* will screen the distance between points when calculating the
    electric field using a gaussian smear.  Default is False (no screening).

    *smear* is the distance (in Bohr) that the grid points are smeared by when
    using screening.  Defaults to 1.0.

    *draw* will draw the electric fields using matplotlib.  If set to False,
    it will instead return the grid points and efields at each grid point.

    '''
    from numpy import array, transpose, arange, zeros, zeros_like, reshape, log10, sqrt, shape, real
    from matplotlib import ticker
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from .constants import BOHR2ANGSTROM, ANGSTROM2BOHR, atomic_radius, atomic_color
    from sys import exit
    #
    # Convert parameters to float arrays
    xpara = array(xpara, dtype='float32')
    ypara = array(ypara, dtype='float32')
    zpara = array(zpara, dtype='float32')
    
    if not "RETARDATION" in chemObj.key:
        exit("needs to be a DIM job run with retardation")
    # Set up screening value
    scrn = 0
    if screen:
        if 'DDA' in chemObj.key:
            scrn = 2
        else:
            scrn = 1

    if not qm:    
        # Convert to Bohr and transpose coordinates
        coords = transpose(chemObj.dim_coordinates*ANGSTROM2BOHR)
        if cg:
            atom_radii = scale*chemObj.dim_cgradii*ANGSTROM2BOHR
        else:
            atom_radii = scale*chemObj.dim_radii()*ANGSTROM2BOHR

        print("Type: ", calctype)
        if dir < 4:
            try:
                print("Direction: ", ['X', 'Y', 'Z'][dir])
            except:
                exit("Must give valid direction index (0, 1, or 2)")

        # for retardation dir should match pol_vec of output file, if not, exit.
        if "RETARDATION" in chemObj.key:
            if dir != chemObj.key["POL_VEC"]:
                exit("Direction supplied must match polarizaton direction specified in calculation, the default pol_vec direction is y")
        else:
            exit("Without retardation effects there is no magnetic field")

        if mag and sq:
            exit("Can only specify E^4, E^2, or E")

        #
        # Collect the coordinates and dipoles depending on the type of DIM
        # calculation chemObj is.
        if calctype == 'ground image':
            dipole = chemObj.dim_dipoles[calctype]
            if chemObj.dim_charges is not None: charges = chemObj.dim_charges[calctype]
        else:
            if dir == 4:
                exit("Need direction of polarization")
            if calctype == 'static scattered':
                dipole = chemObj.dim_dipoles[calctype][dir]
                if chemObj.dim_charges is not None: charges = chemObj.dim_charges[calctype][dir]
            else:
                print("Frequency: ", chemObj.e_frequencies[freq])
                dipole = chemObj.dim_dipoles[calctype][freq][dir]
                if chemObj.dim_charges is not None: charges = chemObj.dim_charges[calctype][freq][dir]

        # if there are no DIM charges (ie a PIM job), the charges will be set to zero
        if chemObj.dim_charges is None: charges = zeros(len(dipole[:,0]))

        # Transpose the dipole array so we can do matrix-vector multiplication
        dipole = transpose(dipole)
        charges = transpose(charges)
        ### QM
    else:
    # With qm==True, calculate fields from hirshfeld partitioned polarizabilities/charges
    # -- Pengchong Oct. 2016
        coords = transpose(chemObj.coordinates*ANGSTROM2BOHR)
        atom_radii=chemObj.radii()
        dip_tot=(chemObj.hirshfeld_induced_dipoles_loc+chemObj.hirshfeld_induced_dipoles_loc)
        dipole=transpose(dip_tot[:,dir,:])    
        charges = zeros(len(dipole[0,:]))

    # Determine which plane we are drawing by checking which parameter is only
    # length 1
    if xpara.size == 1:
        y, z, x, A, B = __build2DGrid(ypara, zpara, xpara)
        title = "X = " + str(xpara)
        labelA = "Y Axis"
        labelB = "Z Axis"
    elif ypara.size == 1:
        x, z, y, A, B = __build2DGrid(xpara, zpara, ypara)
        title = "Y = " + str(ypara)
        labelA = "X Axis"
        labelB = "Z Axis"
    elif zpara.size == 1:
        x, y, z, A, B = __build2DGrid(xpara, ypara, zpara)
        title = "Z = " + str(zpara)
        labelA = "X Axis"
        labelB = "Y Axis"
    else:
        exit("Must define a plane; only parameters can have a range")
    # Initialize the results array to zero and cast is as complex
    if vec:
        # Create a 1D array for each component of the electric field
        # This a the grid flattened
        hfieldx = zeros_like(x)
        hfieldy = zeros_like(x)
        hfieldz = zeros_like(x)
        #tefield = zeros_like(x)
        hfieldx = array(hfieldx, dtype='complex')
        hfieldy = array(hfieldy, dtype='complex')
        hfieldz = array(hfieldz, dtype='complex')
        #tefield = array(efieldz, dtype='complex',order='f')
    else:
        hfield = zeros_like(x)
        hfield = array(hfield, dtype='complex', order='f')
    nsolv = chemObj.key["NSOLV"]
    print('nsolv: ', nsolv)
    if vec:
        from .f2py import drawing_math_vectors_ret_magnetic as dmv_ret_mag
        # Returns efieldx, efieldy, efieldz
        if "RETARDATION" in chemObj.key:
            avec = chemObj.key['A_VEC']
            dmv_ret_mag(x, y, z, coords, dipole, atom_radii, hfieldx, hfieldy, hfieldz,
                       charges, False, dir, scrn, smear, chemObj.e_frequencies[freq], nsolv, avec)
        else:
            exit("you should already be gone")
        # Stack and reshape into grid of vectors
        from numpy import column_stack
        hfield = column_stack((hfieldx, hfieldy, hfieldz))
        hfield = reshape(hfield, (A.shape[0], A.shape[1],3))
        # Return coordinate arrays and the vector array
        return A, B, hfield
    
    from .f2py import drawing_math_ret_magnetic as dm_ret_mag
    if "RETARDATION" in chemObj.key:
        avec = chemObj.key['A_VEC']
        print(chemObj.e_frequencies[freq])
        dm_ret_mag(x, y, z,  coords,  dipole, atom_radii, hfield , charges,
                False, dir, scrn, smear,chemObj.e_frequencies[freq], nsolv, avec) 
    else:
        exit("you should already be gone")
    hfield = reshape(hfield, A.shape)

    if mag:
        hfield = sqrt(hfield)
    elif not sq:
        hfield = hfield*hfield
    #
    # Creates the contour plot
    if draw:
        plt.figure(facecolor='#7F7F7F')
       #plt.axes(axisbg='k')
        pallete = cm.jet
        if max is None:
            cont = plt.contourf(A, B, log10(hfield),
                                locator=ticker.MaxNLocator(500),
                                cmap=pallete,
                                )
        else:
            cont = plt.contourf(A, B, log10(hfield), levels=arange(min,max,step, dtype=float),
                                extend='both',
                                cmap=pallete,
                                )
        CB = plt.colorbar(cont)
        if sq:
            CB.set_label('Magnetic Field Magnitude^2 (log(Atomic Units))')
        else:
            CB.set_label('Magnetic Field Magnitude^4 (log(Atomic Units))')
        plt.title(title + " (Angstroms)")
        plt.xlabel(labelA + " (Angstroms)")
        plt.ylabel(labelB + " (Angstroms)")
        plt.show()
        plt.close()
        del A, B, hfield
        import gc
        gc.collect()
    else:
        return A, B, hfield



def __build2DGrid(aPara, bPara, cPara):
    ''' A generic function that builds a 2 dimensional grid using a and b
        parameters.  It also creates the arrays needed to calculate the
        electric field in space.

    '''
    from .constants import ANGSTROM2BOHR
    from numpy import linspace, reshape, ones_like, meshgrid,round
    from sys import exit
    # Quack
    try:
        minA, maxA, stepA = aPara
        minB, maxB, stepB = bPara
    except TypeError:
        exit(1, "Only on parameter should be an integer!")
    # Determine the number of data points in a single row/column
    tempB = (maxB - minB)/stepB + 1
    tempA = (maxA - minA)/stepA + 1
    # Build the grid
    a = linspace(minA, maxA, int(round(tempA)))
    b = linspace(minB, maxB, int(round(tempB)))
    aa, bb = meshgrid(a,b)
    # Flatten the grid to speed up the looping later
    a = reshape(aa, aa.size)
    b = reshape(bb, bb.size)
    c = ones_like(a)*cPara
    # Convert units
    a = ANGSTROM2BOHR(a)
    b = ANGSTROM2BOHR(b)
    c = ANGSTROM2BOHR(c)
    # Return arrays
    return a, b, c, aa, bb

def raman_draw(cdata, sticks=True, poop=False, width=10.0, scale_freq=1.0,
               dim=(8,6), dpi=300, lw=2.0, fs=20, invert=False, **kwargs):
    from mfunc import sum_lorentzian
    from constants import PI
    import matplotlib.pyplot as plt
    from numpy import linspace
    if 'dir' in kwargs:
        cdata.collect_raman_derivatives(dir=kwargs['dir'], poop=poop)
        cdata._raman = None
    elif 'RAMAN' in cdata.calctype:
        pass
    else:
        cdata.collect_raman_derivatives(poop=poop)

    cdata.scale_vfreq(scale_freq)
    freq  = cdata.v_frequencies
    cross = cdata.cross_section()
    #scale = (2 * 10**36) / (20 * PI)
    print('WARNING: scale change from 10**36 to 10**32')
    scale = (2 * 10**32) / (20 * PI)

    # adding axis Lables
    # the plot is now ax for the subplot
    fig = plt.figure(figsize=dim, dpi=dpi)
    ax = fig.add_subplot(111)
    ax.set_ylabel('$d\sigma/d\Omega\ (10^{-32}cm^2/sr)$', fontsize=fs)
    ax.set_xlabel('wave number (cm$^{-1}$)', fontsize=fs, family='serif')
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    sum = 0
    if sticks:
        print('Frequency     Cross Section')
        for i in xrange(len(cross)):
            print(freq[i], '      ', cross[i]*scale)
            ax.plot((freq[i], freq[i]), (0, cross[i]*scale), 'r')
            sum += cross[i]*scale
        print('Total Cross Sectional Area: ', sum)

    if 'min' in kwargs:
        minx = kwargs['min']
    else:
        minx = 0
    if 'max' in kwargs:
        maxx = kwargs['max']
    else:
        maxx = freq[-1]*1.2

    plt.xlim(minx, maxx)
    domain = linspace(minx, maxx, num=2000)
    # calculate lorentzian spectrum
    y = sum_lorentzian(domain, peak=freq, height=cross, hwhm=width)
    ax.plot(domain, y*(10**32), 'k', lw=lw)
    if invert:
        ax.invert_xaxis()
    if 'file' in kwargs:
        plt.savefig(kwargs['file'], dpi=300, format='png')
    else:
        plt.show()

def vroa_draw(cdata, sticks=True, poop=False, width=10.0, scale_freq=1.0,
               direction='180deg', dim=(4,3), dpi=300, lw=2.0, fs=8, **kwargs):
    from mfunc import sum_lorentzian
    from constants import PI
    import matplotlib.pyplot as plt
    from numpy import linspace
    import numpy
    if 'dir' in kwargs:
        cdata.collect_raman_derivatives(dir=kwargs['dir'], poop=poop)
        cdata._raman = None
    elif 'VROA' in cdata.calctype:
        intensity = numpy.array(cdata.vroa_intensities[direction])
        cdata.vroa_intensities['180deg'] = intensity
        cdata.vroa_cross_section()
        intensity = cdata.vroa_intensities['180deg']
    else:
        cdata.collect_roa_derivatives(poop=poop)
        cdata.calc_roa_intensities()
        intensity = numpy.array(cdata.vroa_intensities[direction])

    cdata.scale_vfreq(scale_freq)
    freq  = cdata.v_frequencies
    # NOTE frequencies are printed x10E3, so divide by 1000
    intensity = intensity / 1000.0
    #scale = (2 * 10**36) / (20 * PI)
#    print('WARNING: scale change from 10**36 to 10**32')
    scale = 1.0

    # adding axis Lables
    # the plot is now ax for the subplot
    fig = plt.figure(figsize=dim, dpi=dpi)
    ax = fig.add_subplot(111)
    ax.set_ylabel('$\Delta d\sigma/d\Omega\ (180^\\circ) [I^R-I^L)] (\\frac{{cm^2}}{{sr}})$', fontsize=fs)
    ax.set_xlabel('wavenumber (cm$^{-1}$)', fontsize=fs, family='serif')
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    sum = 0
    if sticks:
        print('Frequency     Cross Section')
        for i in xrange(len(intensity)):
            print(freq[i], '      ', intensity[i]*scale)
            ax.plot((freq[i], freq[i]), (0, intensity[i]*scale), 'r')
            sum += abs(intensity[i])*scale
        print('Total Cross Sectional Area: ', sum)

    if 'min' in kwargs:
        minx = kwargs['min']
    else:
        minx = 0
    if 'max' in kwargs:
        maxx = kwargs['max']
    else:
        maxx = freq[-1]*1.2

    plt.xlim(minx, maxx)
    domain = linspace(minx, maxx, num=2000)
    # calculate lorentzian spectrum
    y = sum_lorentzian(domain, peak=freq, height=intensity, hwhm=width)
    ax.plot(domain, y*32*PI, 'k', lw=lw)
    if 'file' in kwargs:
        plt.savefig(kwargs['file'], dpi=150, format='png')
    else:
        plt.show()


def write_cube(xyz,xpara,ypara,zpara,density,name):
    from .constants import ANGSTROM2BOHR, BOHR2ANGSTROM
    from numpy import where,zeros
    coords = ANGSTROM2BOHR(xyz.coordinates)
    atmNum=zeros(xyz.natoms)
    index=where(xyz.atoms=="Co")[0]
    atmNum[index]=27
    index=where(xyz.atoms=="Zn")[0]
    atmNum[index]=30
    index=where(xyz.atoms=="O")[0]
    atmNum[index]=8
    index=where(xyz.atoms=="C")[0]
    atmNum[index]=6
    index=where(xyz.atoms=="N")[0]
    atmNum[index]=7
    index=where(xyz.atoms=="H")[0]
    atmNum[index]=1
    index=where(xyz.atoms=="Ag")[0]
    atmNum[index]=47
    index=where(xyz.atoms=="Au")[0]
    atmNum[index]=79

    xstep=(xpara[1]-xpara[0])/(xpara[2]-1)
    ystep=(ypara[1]-ypara[0])/(ypara[2]-1)
    zstep=(zpara[1]-zpara[0])/(zpara[2]-1)

    f = open(name, 'w')
    print('DRAWING CUBE FILE', file=f)
    print('OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z', file=f)
    print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(xyz.natoms, xpara[0], ypara[0], zpara[0]), file=f)
    print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(xpara[2], xstep, 0.0, 0.0), file=f)
    print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(ypara[2], 0.0, ystep, 0.0), file=f)
    print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}'.format(zpara[2], 0.0, 0.0, zstep), file=f)
    for i in xrange(xyz.natoms):
        print('{0:4d}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}'.format(int(atmNum[i]), 0.0, coords[i,0], coords[i,1], coords[i,2]), file=f)

    t = 0
    for i in range(xpara[2]):
        for j in range(ypara[2]):
            for k in range(zpara[2]):
                print('{0:< 13.5E}'.format(density[i,j,k]),end='', file=f)
                if (t % 6 == 5):
                    print('\n', end='', file=f)
                t += 1
    print('\n', end='', file=f)
    f.close()


def calcfield(qmobj, dimobj, index=0, dir=4, scale=0.7978845608028654):
    '''
    This routine is a wrapper routine to calculate the DIM fields at 
    each QM atom and the center of mass of the QM system. 

    For now, the function is only checked to work with 'FD scattered' 
    calctype. Use others at your own risk, or modify accordingly

    *scale* will scale the radius of the DIM atoms to make them their effective
    radius.  The default value is calculated according to
    J. Phys. Chem. C, 2008, 112 (40), 15697; Eq. 15  (0.7978845608028654)
    This value should be used in most cases where screening was used during the DIM calculation.

    Terms from induced charges are not included
    '''
    from .drawing import drawField
    import numpy as np 

    a = qmobj.center_of_mass[0]
    b = qmobj.center_of_mass[1]
    c = qmobj.center_of_mass[2]
    x = [a, a, 1]
    y = [b, b, 1]
    z = c
    fieldarray = np.zeros(len(qmobj.coordinates[:,0]))

    A, B, tmp = drawField(dimobj, x, y, z, calctype='FD scattered', freq=index, dir=dir, sq=True,
                    draw=False, screen=True, smear=3.0, scale=scale)
    efield_com = (tmp.real)

    for i in range(len(fieldarray)):
        a = qmobj.coordinates[i,0]
        b = qmobj.coordinates[i,1]
        c = qmobj.coordinates[i,2]
        x = [a, a, 1]
        y = [b, b, 1]
        z = c
        A, B, tmp = drawField(dimobj, x, y, z, calctype='FD scattered', freq=index, dir=dir, sq=True,
                        draw=False, screen=True, smear=3.0, scale=scale)
        fieldarray[i] = (tmp.real)

    return efield_com, fieldarray
