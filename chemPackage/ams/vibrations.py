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
        elif line == " Index  Atom      ---- Displacements (x/y/z) ----":
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
    if 'NORMAL MODES END'  in indices:
        e = indices['NORMAL MODES END'] - 5
    else:
        self._raise_or_pass('Error locating the end of the normal mode table')
        return

    self.normal_modes = None
    self.v_frequencies = array([], dtype=float)
    # Find where each vibtational frequency table starts
    ar = [i for i, x in enumerate(f[s:e], s) if vib_freq(x)]
    for s in ar:
        # Add vibrational frequencies, however many on this line
        ln = f[s-1].split()[4]
        self.v_frequencies = append(self.v_frequencies, float(ln))
        # Add the normal modes. First time through, define the
        # normal modes array.  Make the normal modes array 3-d
        # (hence the extra brackets) so that append will work
        # properly.
        s += 1
        e = s + self.natoms
        m = array([[x.split()[2:5] for x in f[s:e]]],dtype=float)
        try:
            self.normal_modes = append(self.normal_modes, m, axis=0)
        except ValueError: 
            self.normal_modes = m

    self.nmodes = len(self.v_frequencies)

    # Simultaneously normalize in Bohr then convert to Angstroms
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
    

    
def collect_frequencies_mbh(self, f, indices):
    from ..constants import BOHR2ANGSTROM as B2A
    from numpy.linalg import norm
    import numpy as np

    if 'MBH' in indices:
        s = indices['MBH']
    else:
        self._raise_or_pass('Error locating the normal mode tables')
        return
    if 'MBH END'  in indices:
        e = indices['MBH END']
    else:
        self._raise_or_pass('Error locating the end of the normal mode table')
        return

    self.IR = array([], float)
    self.v_frequencies = array([], dtype=float)
    for i in range(s, e):
        if 'Index' in f[i] and 'Frequency' in f[i] and 'Intensity' in f[i]: ## normal modes, not mass weighted, in Bohr
            ln = f[i].split()[4]
            ir = f[i].strip('\n').split()[7]
            self.v_frequencies = append(self.v_frequencies, float(ln))
            self.IR = append(self.IR, float(ir))
    


    self.nmodes = len(self.v_frequencies)
    self.normal_modes = np.zeros((self.nmodes, self.natoms, 3))
    ii = 0
    for i in range(s, e):
        if 'Index' in f[i] and 'Frequency' in f[i] and 'Intensity' in f[i]: ## normal modes, not mass weighted, in Bohr
            for j in range(self.natoms):
                self.normal_modes[ii][j][0] = float(f[i+1+j].strip('\n').split()[-3])
                self.normal_modes[ii][j][1] = float(f[i+1+j].strip('\n').split()[-2])
                self.normal_modes[ii][j][2] = float(f[i+1+j].strip('\n').split()[-1])
            ii += 1


    for i in range(self.nmodes):
        self.normal_modes[i] *= B2A / norm(self.normal_modes[i].flatten())
 


