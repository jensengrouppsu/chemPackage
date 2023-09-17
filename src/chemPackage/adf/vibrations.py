from __future__ import print_function, division
from numpy import array, append, argsort, row_stack, reshape, sqrt

def collect_frequencies(self, f, indices):
    '''Collect frequencies and IR intensities.'''
    from ..constants import BOHR2ANGSTROM as B2A
    from numpy.linalg import norm

    # Function to detrermine if line is vibrational frequencies
    def vib_freq(line):
        ln = line.split()
        if not ln:
            return False
        elif len(ln) <= 3 and '------------------------' not in ln[0]:
            return True
        else:
            return False

    # The vibrational frequencies and normal modes are placed, in
    # acending order, in tables that contain the movement in the
    # x, y, and z direction for each atom in the system.  There are
    # three tables per page, but sometimes there are less than three
    # tables and this must be accounted for.  
    #
    # The location of the start of each normal mode table is found as
    # where therre are three or less numbers on a line, using the 
    # funtion 'vib_freq' defined at the top of this file. Then the
    # table is defined as ended when a blank line is found.  Each block
    # Is collected within these bounds.
    #
    # The normal modes are collected as a 3-d array, where axis 0 is the
    # different vibrational frequencies, axis 1 is the atoms, and axis 2
    # is the x, y, and z coordinate.
    if 'NORMAL MODES' in indices:
        s = indices['NORMAL MODES']
    else:
        self._raise_or_pass('Error locating the normal mode tables')
        return
    if 'IR' in indices:
        e = indices['IR'] - 13
    else:
        self._raise_or_pass('Error locating the IR intensities')
        return

    self.normal_modes = None
    self.v_frequencies = array([], dtype=float)
    # Find where each vibtational frequency table starts
    ar = [i for i, x in enumerate(f[s:e], s) if vib_freq(x)]
    # Collect from each table
    for s in ar:
        # Add vibrational frequencies, however many on this line
        ln = array(f[s].split(), dtype=float)
        self.v_frequencies = append(self.v_frequencies, ln[:])
        # Add the normal modes. First time through, define the
        # normal modes array.  Make the normal modes array 3-d
        # (hence the extra brackets) so that append will work
        # properly.
        s += 2
        e = s + self.natoms
        m = array([[x.split()[1:4] for x in f[s:e]]],dtype=float)
        try:
            self.normal_modes = append(self.normal_modes, m, axis=0)
        except ValueError: 
            self.normal_modes = m
        # Use try in case there are less than three columns.  
        try:
            m = array([[x.split()[4:7] for x in f[s:e]]], dtype=float)
            self.normal_modes = append(self.normal_modes, m, axis=0)
            try:
                m = array([[x.split()[7:10] for x in f[s:e]]], dtype=float)
                self.normal_modes = append(self.normal_modes, m, axis=0)
            except (ValueError, IndexError):
                pass
        except (ValueError, IndexError):
            pass
    self.nmodes = len(self.v_frequencies)
    # Simultaneously normalize in Bohr then convert to Angstroms
    #for i in xrange(self.nmodes):
    for i in range(self.nmodes):
        self.normal_modes[i] *= B2A / norm(self.normal_modes[i].flatten())

    # The IR intensities are listed multiple times;  we are only
    # interested in the one listed after the modes have been collected.
    s = indices['IR']
    e = self.nmodes + s
    # Zhongwei: In case the frequency is calculated numerically
    try:
       self.IR = array([x.split()[2] for x in f[s:e]], dtype=float)
    except ValueError:
        s = indices['IR'] + 2
        e = self.nmodes + s
        self.IR = array([x.split()[2] for x in f[s:e]], dtype=float)

def collect_raman(self, f, indices):
    '''Collect Raman scattering factors.'''
    from numpy import fastCopyAndTranspose as fcat
    
    # The Raman intensities are listed in symmetry order, so they may
    # not be in numerical order. To account for this, we record the
    # frequencies and intensities, then use the frequencies to sort the
    # intensities after all are collected.
    if 'LIFETIME' in self.subkey:
        if 'FD RAMAN' in indices:
            ar = indices['FD RAMAN']
        else:
            self._raise_or_pass('Error locating FD Raman intensities')
            return
    else:
        if 'STATIC RAMAN' in indices:
            ar = indices['STATIC RAMAN']
        else:
            self._raise_or_pass('Error locating static Raman intensities')
            return
    for s in ar:
        if 'LIFETIME' in self.subkey:
            e = next(i for i,x in enumerate(f[s:], s) if '=========' in x)
            # For complex calculations, the scattering factors are listed
            # in table form.  The total factor is two lines after the
            # frequency.
            tp = []
            for i in range(s, e, 3):
                ln1 = f[i].split()
                ln2 = f[i+2].split()
                tp.append([float(ln1[0]), float(ln2[1])])
            tp = array(tp)
        else:
            e = next(i for i,x in enumerate(f[s:], s) if not x.strip())
            tp = array([x.split() for x in f[s:e]], dtype=float)
            tp = fcat(array([tp[:,0], tp[:,2]]))
        try:
            self._raman = row_stack((self._raman, tp))
        except ValueError:
            self._raman = tp
    index = argsort(self._raman[:,0])
    self._raman = self._raman[index,1]


def collect_vroa(self, f, indices):
    '''Collect the VROA intensities.'''

    from numpy import zeros

    if 'VROA INTENSITIES' in indices and 'FREQUENCIES' in self.calctype:
        s = indices['VROA INTENSITIES']
    elif 'VROA INTENSITIES' in indices:
        return
    else:
        self._raise_or_pass('Error locating VROA intensities')

    lcomplex =  'LIFETIME' in self.subkey
    if lcomplex:
        self.vroa_intensities = {'freq'  : zeros((self.nmodes), dtype=float),
                                 '0deg'  : zeros((self.nmodes), dtype=complex),
                                 '180deg': zeros((self.nmodes), dtype=complex),
                                 'x90deg': zeros((self.nmodes), dtype=complex), 
                                 'z90deg': zeros((self.nmodes), dtype=complex)}
    else:
        self.vroa_intensities = {'freq'  : zeros((self.nmodes), dtype=float),
                                 '0deg'  : zeros((self.nmodes), dtype=float),
                                 '180deg': zeros((self.nmodes), dtype=float),
                                 'x90deg': zeros((self.nmodes), dtype=float),
                                 'z90deg': zeros((self.nmodes), dtype=float)}

    dir = ['0deg', '180deg', 'x90deg', 'z90deg']
    for i in s:
        ipart = 1
        for j in range(self.nmodes*2+3):
            if 'IMAGINARY' in f[i+j]: ipart = 1j
            try:
                line = [float(x) for x in f[i+j].split()[1:6]]
                if len(line) == 5:
                    for x in range(self.nmodes):
                        if round(line[0],1) == round(self.v_frequencies[x],1):
                            mode = x
                    ##### Corrected indent level -- Pengchong Liu, 02/13/2017 #####
                            self.vroa_intensities['freq'][mode]    = line[0]
                            for k in range(len(dir)):
                                self.vroa_intensities[dir[k]][mode] += line[k+1]*ipart
            except ValueError:
                if not lcomplex:
                    break
                elif ipart == 1:
                    continue
                else:
                    break

def collect_dipder(self, f, indices):
       

    #Collect the dipole derivatives
    if 'DIPOLE DERIVATIVES' in indices:
        #Begin the Collection
        ar = indices['DIPOLE DERIVATIVES']
        s = ar + 1
        e = next(i for i, x in enumerate(f[s:], s) if x == '')
        #print([x.split()[1:4] for x in f[s:e]])
        ar = array([x.split()[1:4] for x in f[s:e]], dtype=float)
        self.dgdip = []
        for ln in ar:
            dder = array([ln[0], ln[1], ln[2]])
            self.dgdip = append(self.dgdip, dder)
        #Format dipole dev componets like the output
        self.dgdip = self.dgdip.reshape(e-s, 3)
        calc_dmudq(self, f, indices)

def calc_dmudq(self, f, indices):
        from numpy import hstack,vstack
 
        #Collect the normal modes atomic displacements
        collect_frequencies(self, f, indices)
        normmodes = array([])
        dgtemp = array([])
        temp = 0
        tmparry = array([])
        #Format the normal mode atomic displacements so they are formated like the dipole der
        #Add in the first column so hstack works
        normmodes= append(normmodes, self.normal_modes[0])
        normmodes = vstack(normmodes)
        for x in range(1,self.nmodes):
            tmparry = self.normal_modes[x]
            normmodes = hstack((normmodes, vstack(tmparry.flatten())))
        count = 0
        #Calculated the mass-weighted step size
        sQ = self.step_size(.01)
        
        normmodes = normmodes / 0.529177249 #convert to bohr

        #Convert the dipole derivatives to x, y, and z components
        for x in range(self.nmodes):
            for y in range(3):
                for z in range(self.natoms*3):
                    temp = temp + self.dgdip[z][y]*normmodes[z][x]
                dgtemp = append(dgtemp, temp)
                temp =0
        dgtemp = dgtemp.reshape(self.nmodes, 3)
        #Mass-weigh the dipole der
        for x in range(self.nmodes):
            for y in range(3):          
                dgtemp[x][y] = dgtemp[x][y] * .01 / (sQ[x])
        self.dgdip = dgtemp

def calc_dim_dmudx (self, f, indices):
    import cmath
    self.dgdip = 0
    self.dgdip = []
    if 'DIM DIPOLE IMAG' in indices:
        loop = 2
    else:
        loop = 1
    for a in range(loop):
        #Collect the Real Dipole
        geo = []
        disp = self.natoms* 6 + 1
        if a == 0:
            ar = indices['DIM DIPOLE REAL']
        elif a == 1:
            ar = indices['DIM DIPOLE IMAG']
#    Save! ADF steps in 0.01 in each direction. It may be different
#       in other cases. This will collect the geometry and calc dx
#
#       dx = []
#       geom = []
#       for x in indices['OPTIMIZED GEOMETRY']:
#            for y in range(3):
#                geom = append(geom, f[x+y].split()[5:8])
#        for x in range(0,disp*9-9,19):
#            temp = float(geom[x+9]) - float(geom[x])
#            dx = append(dx, temp)
#
        for i in ar:
            geo  = append(geo, f[i].split()[0:3])
        
        geo = geo.reshape(disp,3)
        dmudx = []
        for x in range(1,disp,2):
            for y in range(3):
                temp = (float(geo[x+1][y]) - float(geo[x][y]))/0.02
                dmudx = append(dmudx, temp)
        dmudx = dmudx.reshape(self.natoms*3 , 3)
        #Convert to bohr
        dmudx = dmudx*0.529177249
        self.dgdip = dmudx
        calc_dmudq(self, f, indices)
        if a ==0:
            real = self.dgdip
        elif a ==1:
            self.dgdip = self.dgdip*1j + real
