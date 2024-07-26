from numpy import array, where, append, arange, zeros
from numpy import row_stack, column_stack, argsort
from chemPackage import collect

def collect_polarizability(self, f, indices):
    '''Drive the collection of polarizabilities from different methods.'''
    __aoresponse(self, f, indices)



def __aoresponse(self, f, indices):
    '''Collect AORESPONSE polarizability.'''

    # This section collections the AORESPONSE routine tensors. Each
    # polarizability is identified by the frequency in Hartrees that
    # it was calculated at.  
    #
    # AORESPONSE can report imaginary polarizability tensors.  If the
    # imaginary tensor is reported, then it is added to the real to
    # make it complex.
    if 'AORESPONSE' in indices:
        ar = indices['AORESPONSE']
        self.e_frequencies = array([], dtype=float)
        for ix in ar:
            # Collect frequency
            ln = f[ix+1].split()
            self.e_frequencies = append(self.e_frequencies, float(ln[2]))
            # Collect tensor.  Different for real and complex.
            if f[ix+14].strip() =="IMAGINARY POLARIZABILITY":
           #if 'LIFETIME' in self.subkey:
                s = ix + 8
                e = ix + 11
                r = array([[x.split() for x in f[s:e]]], dtype=float)
                try:
                    s = ix + 17
                    e = ix + 20
                    i = array([[x.split() for x in f[s:e]]], dtype=float)
                except ValueError:
                    s = ix + 16
                    e = ix + 19
                    i = array([[x.split() for x in f[s:e]]], dtype=float)
                pol = r + i*1j
            else:
                s = ix + 8
                e = ix + 11
                pol = array([[x.split() for x in f[s:e]]], dtype=float)
            # Add this tensor to the array.
            try: 
                self.qm_pol= append(self.qm_pol, pol, axis=0)
            except ValueError:
                self.qm_pol = pol
        self.npol = len(self.e_frequencies)

def collect_opticalrotation(self, f, indices):
    '''Collect the optical rotation data'''
    if 'OPTICAL ROTATION' in indices:
        ar = indices['OPTICAL ROTATION']
        self.e_frequencies = array([], dtype=float)
        for ix in ar:
            # Collect frequency
            ln = f[ix+1].split()
            self.e_frequencies = append(self.e_frequencies, float(ln[2]))
            # Collect tensor.  Different for real and complex.
          # if 'LIFETIME' in self.subkey:
            if f[ix+14] == ' Imaginary Optical rotation tensor:':
                s = ix + 8
                e = ix + 11
                r = array([[x.split() for x in f[s:e]]], dtype=float)
                try:
                    st='        IMAGINARY OPTICAL ROTATION'
                    s = next(i for i, x in enumerate(f[e:], e) if st in x) + 3
                except StopIteration:
                    st=' Imaginary Optical rotation tensor:'
                    s = next(i for i, x in enumerate(f[e:], e) if st in x) + 2
                e = s + 3
                i = array([[x.split() for x in f[s:e]]], dtype=float)
                ord = r + i*1j
            else:
                s = ix + 8
                e = ix + 11
                ord = array([[x.split() for x in f[s:e]]], dtype=float)
            # Generate G-tensor from ord
            Gtensor = - ord * self.e_frequencies[len(self.e_frequencies)-1]
            # Add this tensor to the array.
            try:
                self.ord = append(self.ord, ord, axis=0)
            except ValueError:
                self.ord = ord
            # Add G-tensor to the array
            # NB: This is the G' tensor: G = -iG'
            # which is from the Imaginary part of the dip--mag dip response
            try:
                self.gtensor = append(self.gtensor, Gtensor, axis=0)
            except ValueError:
                self.gtensor = Gtensor
        self.npol = len(self.e_frequencies)
    else:
        self._raise_or_pass('Error locating AORESPONSE optical rotation')


def collect_atensor(self, f, indices):
    '''Collects the dipole-quadrupole tensor.'''

    if 'A-TENSOR' in indices:
        ar = indices['A-TENSOR']
        row = ['x','y','z']
        col = ['xx','xy','xz','yy','yz','zz']
        for ix in ar:
            imag = False
        #   if 'LIFETIME' in self.subkey:
            if f[ix+19] == ' Imaginary DIPOLE-QUADRUPOLE Polarizability tensor:':
                temp = zeros((1,3,6),dtype=complex)
            else:
                temp = zeros((1,3,6),dtype=float)
            for i in range(39):
                s = ix + i
                string = f[s]
                # cycle if there isn't something to collect
                if '---' in string: continue
                if len(string.split()) <= 1: continue
                # check if we've reached the end of the A-tensors
            #   if 'LIFETIME' in self.subkey:
                if f[ix+19] == ' Imaginary DIPOLE-QUADRUPOLE Polarizability tensor:':
                    if len(string.split()) > 3 and imag: break
                else:
                    if len(string.split()) > 3: break
                # check if we're collecting the real of imaginary component
                if '      IMAGINARY DIPOLE-QUADRUPOLE POLARIZABILITY' in string \
                or ' Imaginary DIPOLE-QUADRUPOLE Polarizability tensor' in string:
                    imag = True
                    continue
                # replace some characters (which can differ)
                remove = [',', 'A', '_', '=']
                for char in remove:
                    string = string.replace(char, ' ')
                # split string and collect the tensor value
                string = string.split()
                try:
                    tempval = float(string[2])
                except ValueError:
                    continue
                except IndexError:
                    continue
                # figure out row and column from string[0] and [1]
                r = [r for r in range(3) if row[r] in string[0]]
                c = [c for c in range(6) if col[c] in string[1]]
                if not imag:
                    temp[0,r,c] += tempval
                else:
                    temp[0,r,c] += tempval*1j
            # add this tensor to the array
            try:
                self.atensor = append(self.atensor, temp, axis=0)
            except ValueError:
                self.atensor = temp
    else:
        self._raise_or_pass('Error locating AORESPONSE DIPOLE-QUADRUPOLE tensor')

def collect_ctensor(self, f, indices):
    '''Collects the quadrupole-quadrupole tensor.'''

    if 'C-TENSOR' in indices:
        ar = indices['C-TENSOR']
        row = ['xx','xy','xz','yy','yz','zz']
        col = ['xx','xy','xz','yy','yz','zz']
        self.ctensor is None
        for ix in ar:
            imag = False
        #   if 'LIFETIME' in self.subkey:
            if f[ix+37] == ' Imaginary QUADRUPOLE-QUADRUPOLE Polarizability tensor:':
                temp = zeros((1,6,6),dtype=complex)
            else:
                temp = zeros((1,6,6),dtype=float)
            for i in range(75):
                s = ix + i
                string = f[s]
                # cycle if there isn't something to collect
                if '---' in string: continue
                if len(string.split()) <= 1: continue
                # check if we've reached the end of the A-tensors
                if f[ix+37] == ' Imaginary QUADRUPOLE-QUADRUPOLE Polarizability tensor:':
                    if len(string.split()) > 3 and imag: break
                else:
                    if len(string.split()) > 3: break
                # check if we're collecting the real of imaginary component
                if ' Imaginary QUADRUPOLE-QUADRUPOLE Polarizability tensor:' in string:
                    imag = True
                    continue
                # replace some characters (which can differ)
                remove = [',', 'A', '_', '=']
                for char in remove:
                    string = string.replace(char, ' ')
                # split string and collect the tensor value
                string = string.split()
                try:
                    tempval = float(string[2])
                except ValueError:
                    continue
                except IndexError:
                    continue
                # figure out row and column from string[0] and [1]
                r = [r for r in range(6) if row[r] in string[0]]
                c = [c for c in range(6) if col[c] in string[1]]
                if not imag:
                    temp[0,r,c] += tempval
                else:
                    temp[0,r,c] += tempval*1j

            # add this tensor to the array
            if self.ctensor is None:
                self.ctensor = array(temp)
            else:
                self.ctensor = append(self.ctensor, temp, axis=0)
    else:
        self._raise_or_pass('Error locating AORESPONSE QUADRUPOLE-QUADRUPOLE tensor')


def collect_magnetizability(self, f, indices):
    '''Collects the magnetic dipole-magnetic dipole polarizability.'''

    if 'MAGNETIZABILITY' in indices:
        ar = indices['MAGNETIZABILITY']
        self.magnetizability = None
        for ix in ar:
            imag = False
           #if 'LIFETIME' in self.subkey:
            if f[ix+14] == ' Imaginary Paramagnetic magnetizability tensor:':
                temp = zeros((1,3,3),dtype=complex)
                vals = [0,14]
            else:
                temp = zeros((1,3,3),dtype=float)
                vals = [0]
            for i in vals:
                s = ix + i
                string = f[s]
                # check if we're collecting the real of imaginary component
                if ' Imaginary Paramagnetic magnetizability tensor:' in string:
                    imag = True
                for j in range(3):
                    string = f[s+2+j]
                    string = string.split()
                    for k in range(3):
                        string[k] = float(string[k])
                        if imag: string[k] = string[k] * 1.j
                        temp[0][j][k] += string[k]
            # add this tensor to the array
            if self.magnetizability is None:
                self.magnetizability = array(temp)
            else:
                self.magnetizability = append(self.magnetizability, temp, axis=0)
    else:
        self._raise_or_pass('Error locating AORESPONSE MAGNETIC DIPOLE-MAGNETIC DIPOLE tensor')



def collect_linearresponse(self, f, indices):
    ''' collects the condensed linear repsonse function'''
    
    if 'LINEAR RESPONSE' in indices:
        temp = zeros((1,self.natoms, self.natoms), dtype=float)
        ar = indices['LINEAR RESPONSE']
        for ix in ar:
            i = 0
            atomcounter = -1
            while len(f[ix+i]) != 0:
                string = f[ix+i]
                stringS = string.split()
                if ':' in string:
                    atomcounter = atomcounter +1
                    for j in range(len(stringS)-3):
                        temp[0,j,atomcounter] = stringS[j+3]
                        temp[0,atomcounter,j] = stringS[j+3]
                    lincontiunation = 10
                else :
                    for j in range(len(stringS)):
                        temp[0,j+lincontiunation,atomcounter] = stringS[j]
                        temp[0,atomcounter,j+lincontiunation] = stringS[j]
                    lincontiunation = lincontiunation +10
                    
                i = i +1
        

            try:
                self.linresp = append(self.linresp, temp, axis=0)
            except ValueError:
                self.linresp = temp
        
