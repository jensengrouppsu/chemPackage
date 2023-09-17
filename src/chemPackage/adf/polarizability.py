from __future__ import print_function, division
from numpy import array, where, append, arange, zeros
from numpy import row_stack, column_stack, argsort
from chemPackage import collect

def collect_polarizability(self, f, indices):
    '''Drive the collection of polarizabilities from different methods.'''
    if 'RAMAN' in self.subkey:
        __polarizability_derivatives(self, f, indices)
        __tensor_derivatives(self, f, indices)
    elif 'RESPONSE' in self.key:
        __response(self, f, indices)
    else:
        __aoresponse(self, f, indices)
        #__atomic(self, fo)

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
            if 'LIFETIME' in self.subkey:
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
            if 'LIFETIME' in self.subkey:
                temp = zeros((1,3,6),dtype=complex)
            else:
                temp = zeros((1,3,6),dtype=float)
            for i in range(42):
                s = ix + i
                string = f[s]
                # cycle if there isn't something to collect
                if '---' in string: continue
                if len(string.split()) <= 1: continue
                # check if we've reached the end of the A-tensors
                if 'LIFETIME' in self.subkey:
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
            if 'LIFETIME' in self.subkey:
                temp = zeros((1,6,6),dtype=complex)
            else:
                temp = zeros((1,6,6),dtype=float)
            for i in range(80):
                s = ix + i
                string = f[s]
                # cycle if there isn't something to collect
                if '---' in string: continue
                if len(string.split()) <= 1: continue
                # check if we've reached the end of the A-tensors
                if 'LIFETIME' in self.subkey:
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


def collect_btensor(self, f, indices):
    '''Collects the dipole-octupole polarizabilities.'''

    if 'B-TENSOR' in indices:
        ar = indices['B-TENSOR']
        row = ['x','y','z']
        col = ['xxx','xxy','xxz','xyy','xyz','xzz','yyy','yyz','yzz','zzz']
        self.btensor = None
        for ix in ar:
            imag = False
            if 'LIFETIME' in self.subkey:
                temp = zeros((1,3,10),dtype=complex)
            else:
                temp = zeros((1,3,10),dtype=float)
            for i in range(64):
                s = ix + i
                string = f[s]
                # cycle if there isn't something to collect
                if len(string.split()) <= 1: continue
                # check if we've reached the end of the A-tensors
                if 'LIFETIME' in self.subkey:
                    if len(string.split()) > 3 and imag: break
                else:
                    if len(string.split()) > 3: break
                # check if we're collecting the real of imaginary component
                if ' Imaginary DIPOLE-OCTUPOLE Polarizability tensor:' in string:
                    imag = True
                    continue
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
                c = [c for c in range(10) if col[c] in string[1]]
                if not imag:
                    temp[0,r,c] += tempval
                else:
                    temp[0,r,c] += tempval*1j
            # add this tensor to the array
            if self.btensor is None:
                self.btensor = array(temp)
            else:
                self.btensor = append(self.btensor, temp, axis=0)
    else:
        self._raise_or_pass('Error locating AORESPONSE DIPOLE-OCTUPOLE tensor')


def collect_dtensor(self, f, indices):
    '''Collects the quadruopole-magnetic dipole polarizability.'''

    if 'D-TENSOR' in indices:
        ar = indices['D-TENSOR']
        row = ['xx','xy','xz','yy','yz','zz']
        col = ['x','y','z']
        self.dtensor = None
        for ix in ar:
            imag = False
            if 'LIFETIME' in self.subkey:
                temp = zeros((1,6,3),dtype=complex)
            else:
                temp = zeros((1,6,3),dtype=float) 
            for i in range(40):
                s = ix + i
                string = f[s]
                # cycle if there isn't something to collect
                if len(string.split()) <= 1: continue
                # check if we've reached the end of the D-tensor
                if 'LIFETIME' in self.subkey:
                  # if len(string.split()) <> 3 and imag: break
                    if len(string.split()) <= 3 and imag: break
                else:
                  # if len(string.split()) <> 3: break
                    if len(string.split()) <= 3: break
                # check if we're collecting the real of imaginary component
                if ' Imaginary QUADRUPOLE-MAGNETIC DIPOLE Polarizability tensor:' in string:
                    imag = True
                    continue
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
                c = [c for c in range(3) if col[c] in string[1]]
                if not imag:
                    temp[0,r,c] += tempval
                else:
                    temp[0,r,c] += tempval*1j
            # add this tensor to the array
            if self.dtensor is None:
                self.dtensor = array(temp)
            else:
                self.dtensor = append(self.dtensor, temp, axis=0)
        # Multiply by 3/2 since the traceless D-tensor are NOT printed in ADF
        self.dtensor = self.dtensor * 1.5
    else:
        self._raise_or_pass('Error locating AORESPONSE QUADRUPOLE-MAGNETIC DIPOLE tensor')


def collect_magnetizability(self, f, indices):
    '''Collects the magnetic dipole-magnetic dipole polarizability.'''

    if 'MAGNETIZABILITY' in indices:
        ar = indices['MAGNETIZABILITY']
        self.magnetizability = None
        for ix in ar:
            imag = False
            if 'LIFETIME' in self.subkey:
                temp = zeros((1,3,3),dtype=complex)
                vals = [0]
            else:
                temp = zeros((1,3,3),dtype=float)
                vals = [0,14]
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


def __response(self, f, indices):
    '''Collect RESPONSE polarizability.'''

    # This section collects data from the RESPONSE routine.  Each
    # polarizability is identified by the frequency in Hartrees that
    # it was calculated at.  The response tensor is listed as YZX,
    # so it has to be reordered to be stored properly.
    if 'RESPONSE' in indices:
        s = indices['RESPONSE']
        e = indices['ADF EXIT'] if 'ADF EXIT' in indices else len(f)
        ar = [i for i, x in enumerate(f[s:e], s) if 'FREQUENCY CYCLE NR.' in x]
        self.e_frequencies = array([], dtype=float)
        # Index to reorder tenosr
        i = array([2, 0, 1])
        for ix in ar:
            # Collect Frequency
            ln = f[ix].split()
            self.e_frequencies = append(self.e_frequencies, float(ln[5]))
            # Collect tensor
            s = ix + 2
            e = ix + 5
            pol = array([x.split() for x in f[s:e]], dtype=float)
            # RESPONSE tensor is YZX, so reorganize and add
            try:
                self.qm_pol = append(self.qm_pol, array([pol[i][:,i]]), axis=0)
            except ValueError:
                self.qm_pol = array([pol[i][:,i]])
        self.npol = len(self.e_frequencies)
    else:
        self._raise_or_pass('Error locating RESPONSE header')


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
            if 'LIFETIME' in self.subkey:
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
    else:
        if array(['EL_DIPOLE_EL_DIPOLE' in i.upper() for i in self.key['AORESPONSE']]).any():
            self._raise_or_pass('Error locating AORESPONSE polarizabilities')
        else:
            # The dipole-dipole polarizability may not be calculated if we use the new
            # pol code.
            self.calctype.remove('POLARIZABILITY')


def __polarizability_derivatives(self, f, indices):
    'Collect polarizability derivatives from a RAMAN calculation.'

    if 'POLARIZABILITY DERIVATIVES' in indices:
        ar = indices['POLARIZABILITY DERIVATIVES']
        for ix in ar:
            # Collect derivatives differently for real and complex
            if 'LIFETIME' in self.subkey:
                # Find the start end ending index of these derivatives
                # Real part
                s = next(i for i, x in enumerate(f[ix:], ix) if 'Real' in x)
                s = s + 4
                e = next(i for i, x in enumerate(f[s:], s) if not x.strip())
                ar = array([x.split() for x in f[s:e]], dtype=float)
                # Imaginary part
                s = next(i for i, x in enumerate(f[e:],e) if 'Imaginary' in x)
                s = s + 4
                e = next(i for i, x in enumerate(f[s:], s) if '========' in x)
                im = array([x.split() for x in f[s:e]], dtype=float)
                ar = ar + im*1j
            else:
                s = ix + 5
                e = next(i for i, x in enumerate(f[s:], s) if '========' in x)
                ar = array([x.split() for x in f[s:e]], dtype=float)
            # Collect the polarizability derivatives for these normal modes
            for ln in ar:
                pol = array([[[ln[1], ln[2], ln[5]],
                              [ln[2], ln[3], ln[4]],
                              [ln[5], ln[4], ln[6]]]])
                try:
                    self.qm_pol= append(self.qm_pol, pol, axis=0)
                except ValueError:
                    self.qm_pol = pol
            # Collect vibrational frequency
            try:
                nmodes = append(nmodes, ar[:,0].real)
            except UnboundLocalError:
                nmodes = ar[:,0].real
        index = argsort(nmodes)
        self.qm_pol = self.qm_pol[index]
        self.npol = len(self.v_frequencies)

        # set incident frequencies of the calculation
        self.e_frequencies = zeros(self.npol) 
        # keep the zeros for static calculation
        # if freq. dependent calculation collect frequencies
        #   we will collect first freq. and assume 
        #   all other freqs are the same this should be 
        #   the case for a RAMAN calculation
        if 'AORESPONSE' in indices:
            if 'LIFETIME' in self.subkey: 
                i  = indices['AORESPONSE']
                ln = f[i[0]+1].split()               
                self.e_frequencies = self.e_frequencies + float(ln[2])
        elif 'RESPONSE' in indices:
            i  = indices['RESPONSE']
            ln = f[i+2].split()                                                 
            self.e_frequencies = self.e_frequencies + float(ln[5]) 


    else:
        self._raise_or_pass('Error locating polarizability derivatives')


def __atomic(self, fo):
    '''Collect the atomic polarizabilities.'''

    # These actually appear in the output file before the frequency and
    # total tensor is given. The atomic polarizabilites are listed as all
    # the x-row, y-row, and z-row values separately, so we must form the
    # tensors on the fly.
    #
    # The location of each x, y, and z block is located before hand, and
    # these are traversed in order so that we can collect the atomic
    # tensors by frequency so that they line up with the total
    # polarizability tensors.  A clever use of numpy functions 'zips' the
    # blocks for each frequency together into the 3x3 atomic tensors.
    sl = ['     ATOM      xx         xy         xz',
          '     ATOM      yx         yy         yz',
          '     ATOM      zx         zy         zz']
    try:
        arx = where(fo == sl[0])[0]
        ary = where(fo == sl[1])[0]
        arz = where(fo == sl[2])[0]
        test = arz[0]
    except IndexError:
        pass
    else:
        # Create generators to locate each atomic polarizability block
        fnx = (i for i in arx)
        fny = (i for i in ary)
        fnz = (i for i in arz)
        ix = next(fnx, None) # Grab x position to start
        while ix: # Termination is when there are no more x tables
            # Collect atomic polarizabilities for this frequency
            s = ix + 1
            e = next(i for i, x in enumerate(fo[s:], s) if '----' in x)
            xp = [x.split()[2:5] for x in fo[s:e]]
            s = next(fny) + 1
            e = next(i for i, y in enumerate(fo[s:], s) if '----' in y)
            yp = [y.split()[2:5] for y in fo[s:e]]
            s = next(fnz) + 1
            e = next(i for i, z in enumerate(fo[s:], s) if '----' in z)
            zp = [z.split()[2:5] for z in fo[s:e]]
            for i in range(len(xp)):
                for j in range(3):
                    if '*****' in xp[i][j]: xp[i][j] = float('inf')
                    if '*****' in yp[i][j]: yp[i][j] = float('inf')
                    if '*****' in zp[i][j]: zp[i][j] = float('inf')
            xp = array(xp, dtype=float)
            yp = array(yp, dtype=float)
            zp = array(zp, dtype=float)
            # Place the polarizabilities into tensors.
            # The below set of functions will take the three stacks and
            # organize them into tensors by atom.
            atmpol = column_stack((xp,yp,zp))
            try:
                self.polarizability_atomic = row_stack((
                                 self.polarizability_atomic, atmpol))
            except ValueError:
                self.polarizability_atomic = atmpol
            # Find next x table
            ix = next(fnx, None)
        self.polarizability_atomic = (
                self.polarizability_atomic.reshape(self.npol, self.natoms,3,3))


def collect_hyperpolarizability(self, f, indices):
    '''Collect hyperpolarizability.'''

    # This line should appear twice, find both
    if 'HYPERPOLARIZABILITY' in indices:
        indx = indices['HYPERPOLARIZABILITY']
        self.hyperpolarizability = {}

        # Enumerate shortcut
        en = enumerate
        # Locate further hyperpolarizability headers, if any
        # Only search between two instances of the sl line
        word = 'hyperpolarizability tensor beta'
        ar = [i for i, x in en(f[indx[0]:indx[1]], indx[0]) if word in x]
        # There should at least be the sl line in the above search
        for i, ix in en(ar):

            # Find the key.  It is the word in caps
            key = f[ix].split()[1]           
 
            # First index is where 'beta(' first appears.
            try:
                s = next(i for i, x in en(f[ix:ar[i+1]], ix) if 'beta(' in x)
            except IndexError:
                try:
                    s = next(i for i, x in en(f[ix:], ix) if 'beta(' in x)
                except StopIteration:
                    # This occurs if the 'next' expression finds nothing, 
                    # and means all values are zero, since the program won't
                    # print a zero value (which is a stupid programming
                    # decision)
                    self.hyperpolarizability[key] = zeros((3,3,3))
                    continue
            except StopIteration:
                self.hyperpolarizability[key] = zeros((3,3,3))
                continue
            

            # Correlation dictionary
            corr = { 'x' : 0, 'y' : 1, 'z' : 2 }
            # Find last beta value
            e = next(i for i, x in en(f[s:], s) if x[0:5] != 'beta(')
            # Split each line on the '=', then keep 1st number on right of
            # the '='.  Determine where to place the number by the
            # directions listed on the left of '='.
            self.hyperpolarizability[key] = zeros((3,3,3))
            for ln in f[s:e]:
                d, hp = ln.split('=')
                hp = float(hp.split()[0])
                # Extract the three directions as indices
                ix = [corr[d[5]], corr[d[7]], corr[d[9]]]
                # Place into array at correct location
                self.hyperpolarizability[key][ix[0], ix[1], ix[2]] = hp
    else:
        self._raise_or_pass('Error locating hyperpolarizability tensors')

    # Time to grab the invarients
    if 'HYPERPOLARIZABILITY INVARIENTS' in indices:
        end = indices['HYPERPOLARIZABILITY INVARIENTS']
        self.hyperinvarients = {}

        # Locate further hyperpolarizability headers, if any
        # Only search between two instances of the sl line
        word = 'hyperpolarizability tensor beta'
        ar = [i for i, x in en(f[indx[1]:end], indx[1]) if word in x]
        # There should at least be the sl line in the above search
        for i, ix in en(ar):

            # Find the key.  It is the word in caps
            key = f[ix].split()[1]
            # Initiallize the dict
            self.hyperinvarients[key] = { 'beta'     : zeros(3, dtype=float),
                                          'mu'       : zeros(3, dtype=float),
                                          'beta_bar' : 0.0,
                                          'beta_vec' : 0.0,
                                        }
 
            # Grab the "beta" and "mu" variables.  Some of these may not
            # be printed in every calculation.
            # z-component
            self.hyperinvarients[key]['beta'][2] = float(f[ix+1].split()[2])
            self.hyperinvarients[key]['mu'][2]   = float(f[ix+1].split()[5])
            # y-component
            if 'beta_y' in f[ix+2]:
                self.hyperinvarients[key]['beta'][1] = float(f[ix+2].split()[2])
                self.hyperinvarients[key]['mu'][1]   = float(f[ix+2].split()[5])
            else:
                self.hyperinvarients[key]['beta'][1] = float(0.00)
                self.hyperinvarients[key]['mu'][1]   = float(0.00)
            # x-component
            if 'beta_x' in f[ix+3]:
                self.hyperinvarients[key]['beta'][0] = float(f[ix+3].split()[2])
                self.hyperinvarients[key]['mu'][0]   = float(f[ix+3].split()[5])
            else:
                self.hyperinvarients[key]['beta'][0] = float(0.00)
                self.hyperinvarients[key]['mu'][0]   = float(0.00)

            # Grab "beta_bar" anf "beta_vec"
            if 'beta_bar' in f[ix+4]:
                self.hyperinvarients[key]['beta_bar'] = float(f[ix+4].split()[3])
                self.hyperinvarients[key]['beta_vec'] = float(f[ix+4].split()[4])
            elif 'beta_bar' in f[ix+3]:
                self.hyperinvarients[key]['beta_bar'] = float(f[ix+3].split()[3])
                self.hyperinvarients[key]['beta_vec'] = float(f[ix+3].split()[4])
            else:
                self.hyperinvarients[key]['beta_bar'] = float(f[ix+2].split()[3])
                self.hyperinvarients[key]['beta_vec'] = float(f[ix+2].split()[4])
    else:
        self._raise_or_pass('Error locating hyperpolarizability invarients')

# Collect the beta from aoresponse calculations(In this case, from 2n+1 specifically)
def collect_aoresponse_hyperpolarizability(self, f, indices):
    '''Collect aoresponse hyperpolarizability.'''

    # See read_file.py, set up a new term in dictionary
    if 'HYPERPOLARIZABILITY' in indices:
        # Set up the range for scanning, s: start; e: end.
        # Start right at the first line,use [0] due to numerical calculation(e = s + 27)
        s = indices['HYPERPOLARIZABILITY'][0]
        e = s + 27
        self.hyperpolarizability = {}

        # Make an empty list for the real part of beta
        r = []
        # Make an empty list for the imaginary part of beta 
        if 'LIFETIME' in self.subkey:
            i = []
        # Scan the lines from s(tart) to e(nd) by using a loop  
        for j in range(s, e):
            # Return a list of the words in the string, name it as 'o'
            o = f[j].split()
            # Add the fourth element of 'o' to 'r'(***Python starts at 0)
            # Zhongwei: in case there are too many digits in a number
            if '-' in o[3][1:-1]:
               imag = o[3][15:-1]
               o[3] = o[3][0:15]
            r.append(o[3])
            # Add the fifth element of 'o' to 'i'
            if 'LIFETIME' in self.subkey:
               try:
                  i.append(o[4])
               except IndexError:
                  i.append(imag)

        if 'LIFETIME' in self.subkey:
            # Split the line that starts with omega1,omega2 in the output
            p = f[s-4].split()
            # Make an array for the two frequencies
            self.b_e_frequencies = array([], dtype=float)
            self.c_e_frequencies = array([], dtype=float)
            # Fill the list/array with the cooresponding frequency
            self.b_e_frequencies = append(self.b_e_frequencies, float(p[1]))
            self.c_e_frequencies = append(self.c_e_frequencies, float(p[2]))

            # Make an array for 'r', 'i' and reshape each of them
            rhpoltmp = array(r, dtype=float)
            rhpoltmp = rhpoltmp.reshape(3,3,3)
            ihpoltmp = array(i, dtype=float)
            ihpoltmp = ihpoltmp.reshape(3,3,3)
            # Combine real and imaginary part together
            hpoltmp = rhpoltmp + ihpoltmp*1j
            self.hyperpolarizability['FD'] = hpoltmp
        else:
                q = f[s-4].split()
                # Convert string to number
                freq1 = float(q[1])
                freq2 = float(q[2])
                if freq1 == 0 and freq2 == 0:
                    hpoltmp = array(r, dtype=float)
                    hpoltmp = hpoltmp.reshape(3,3,3)
                    self.hyperpolarizability['STATIC'] = hpoltmp
                elif freq1 == -freq2:
                    hpoltmp = array(r, dtype=float)
                    hpoltmp = hpoltmp.reshape(3,3,3)
                    self.hyperpolarizability['OR'] = hpoltmp
                elif freq1 == freq2:
                    hpoltmp = array(r, dtype=float)
                    hpoltmp = hpoltmp.reshape(3,3,3)
                    self.hyperpolarizability['SHG'] = hpoltmp
                else:
                    hpoltmp = array(r, dtype=float)
                    hpoltmp = hpoltmp.reshape(3,3,3)
                    self.hyperpolarizability['EOPE'] = hpoltmp

# Collect the gamma from aoresponse calculations(In this case, from 2n+1 specifically)
def collect_aoresponse_second_hyperpolarizability(self, f, indices):
    '''Collect aoresponse second hyperpolarizability.'''

    # See read_file.py, set up a new term in dictionary
    if 'SECOND HYPERPOLARIZABILITY' in indices:
        # Set up the range for scanning, s: start; e: end.
        # Start right at the first line,use [0] due to numerical calculation(e = s + 81)
        s = indices['SECOND HYPERPOLARIZABILITY'][0]
        e = s + 81
        self.secondhyperpolarizability = {}

        # Make an empty list for the real part of beta
        r = []
        # Make an empty list for the imaginary part of beta 
        if 'LIFETIME' in self.subkey:
            i = []
        # Scan the lines from s(tart) to e(nd) by using a loop  
        for j in range(s, e):
            # Return a list of the words in the string, name it as 'o'
            o = f[j].split()
            # Add the fifth element of 'o' to 'r'(***Python starts at 0)
            r.append(o[4])
            # Add the sixth element of 'o' to 'i'
            if 'LIFETIME' in self.subkey:
                i.append(o[5])

        if 'LIFETIME' in self.subkey:
            # Split the line that starts with omega1,omega2,omega3 in the output
            p = f[s-4].split()
            # Make an array for the three frequencies
            self.b_e_frequencies = array([], dtype=float)
            self.c_e_frequencies = array([], dtype=float)
            self.d_e_frequencies = array([], dtype=float)
            # Fill the list/array with the cooresponding frequency
            self.b_e_frequencies = append(self.b_e_frequencies, float(p[1]))
            self.c_e_frequencies = append(self.c_e_frequencies, float(p[2]))
            self.d_e_frequencies = append(self.d_e_frequencies, float(p[3]))

            # Make an array for 'r', 'i' and reshape each of them
            rshpoltmp = array(r, dtype=float)
            rshpoltmp = rshpoltmp.reshape(3,3,3,3)
            ishpoltmp = array(i, dtype=float)
            ishpoltmp = ishpoltmp.reshape(3,3,3,3)
            # Combine real and imaginary part together
            shpoltmp = rshpoltmp + ishpoltmp*1j
            self.secondhyperpolarizability['FD'] = shpoltmp
        else:
                q = f[s-4].split()
                # Convert string to number
                freq1 = float(q[1])
                freq2 = float(q[2])
                freq3 = float(q[3])
                if freq1 == 0 and freq2 == 0 and freq3 == 0:
                    shpoltmp = array(r, dtype=float)
                    shpoltmp = shpoltmp.reshape(3,3,3,3)
                    self.secondhyperpolarizability['STATIC'] = shpoltmp
                elif freq1 == freq2 == freq3:
                    shpoltmp = array(r, dtype=float)
                    shpoltmp = shpoltmp.reshape(3,3,3,3)
                    self.secondhyperpolarizability['THG'] = shpoltmp
                elif freq1 == freq2 == -freq3:
                    shpoltmp = array(r, dtype=float)
                    shpoltmp = shpoltmp.reshape(3,3,3,3)
                    self.secondhyperpolarizability['IDRI'] = shpoltmp
                elif freq1 == freq2 and freq3 == 0 :
                    shpoltmp = array(r, dtype=float)
                    shpoltmp = shpoltmp.reshape(3,3,3,3)
                    self.secondhyperpolarizability['EFISHG'] = shpoltmp
                elif freq1 == 0 and freq2 == -freq3:
                    shpoltmp = array(r, dtype=float)
                    shpoltmp = shpoltmp.reshape(3,3,3,3)
                    self.secondhyperpolarizability['EFIOR'] = shpoltmp
                else:
                    shpoltmp = array(r, dtype=float)
                    shpoltmp = shpoltmp.reshape(3,3,3,3)
                    self.secondhyperpolarizability['OKE'] = shpoltmp

def __tensor_derivatives(self, f, indices):
    '''Collects all the tensor derivatives, including the polarizability
    derivatives printed in a different format.'''

    scomplex = 'float'
    if 'LIFETIME' in self.subkey: scomplex = 'complex'

    def grab_derivatives(tensor, ix, factor=1.0, dim=3):
        '''Since all the derivatives are printed in a similar fashion, we
        can use a very general algorithm to collect them all.'''

        def tindex(string):
            '''Returns a key index based on a string.'''
            if len(string) == 1:
                keys = ['X','Y','Z']
                return keys.index(string)
            else:
                keys = ['X,X','X,Y','X,Z','Y,X','Y,Y','Y,Z','Z,X','Z,Y','Z,Z']
                a = keys.index(string)
                return a/3, a%3

        for i in range(dim*len(self.v_frequencies)):
            line = f[ix+i]
            if '=====' in line: break
            line = line.split()
            if len(line) == 5:
                key = '{0:0.2f}'.format(float(line[0]))
                del line[0]
            if dim == 3:
                j = tindex(line[0])
                for k in range(3):
                    tensor[key][j][k] += float(line[k+1]) * factor
            else:
                j, k = tindex(line[0])
                for l in range(3):
                    tensor[key][j][k][l] += float(line[l+1]) * factor
        return tensor
            

    if 'POLARIZABILITY DERIVATIVES R' in indices:
        # initialize the polarizability derivative attribute
        self.polarizability_derivatives = {}
        for freq in self.v_frequencies:
            self.polarizability_derivatives['{0:0.2f}'.format(freq)] = zeros((3,3),dtype=scomplex)
        ar = indices['POLARIZABILITY DERIVATIVES R']
        for ix in ar:
            self.polarizability_derivatives = grab_derivatives(self.polarizability_derivatives,
                                                               ix)
        if 'POLARIZABILITY DERIVATIVES I' in indices:
            ar = indices['POLARIZABILITY DERIVATIVES I']
            for ix in ar:
                self.polarizability_derivatives = grab_derivatives(self.polarizability_derivatives,
                                                               ix, factor=1.j)
        
    if 'OPTICAL ROTATION DERIVATIVES R' in indices:
        # initialize the optical rotation derivative attribute
        self.opt_rot_derivatives = {}
        for freq in self.v_frequencies:
            self.opt_rot_derivatives['{0:0.2f}'.format(freq)] = zeros((3,3),dtype=scomplex)
        ar = indices['OPTICAL ROTATION DERIVATIVES R']
        for ix in ar:
            self.opt_rot_derivatives = grab_derivatives(self.opt_rot_derivatives, ix)
        if 'OPTICAL ROTATION DERIVATIVES I' in indices:
            ar = indices['OPTICAL ROTATION DERIVATIVES I']
            for ix in ar:
                self.opt_rot_derivatives = grab_derivatives(self.opt_rot_derivatives,
                                                            ix, factor=1.j)

    if 'A-TENSOR DERIVATIVES R' in indices:
        # initialize the A-tensor derivative attribute
        self.atensor_derivatives = {}
        for freq in self.v_frequencies:
            self.atensor_derivatives['{0:0.2f}'.format(freq)] = zeros((3,3,3),dtype=scomplex)
        ar = indices['A-TENSOR DERIVATIVES R']
        for ix in ar:
            self.atensor_derivatives = grab_derivatives(self.atensor_derivatives, ix, dim=9)
        if 'A-TENSOR DERIVATIVES I' in indices:
            ar = indices['A-TENSOR DERIVATIVES I']
            for ix in ar:
                self.atensor_derivatives = grab_derivatives(self.atensor_derivatives,
                                                            ix, factor=1.j, dim=9)

def collect_quadrupole_beta(self, f, indices):
    '''Collects all the quadrupole beta tensors.'''

    lcmplx = [True if ('FD' in self.calctype or 'LIFETIME' in self.subkey) else False][0]
    dtype = [complex if lcmplx else float][0]
    ip = [2 if lcmplx else 0][0]

    if 'HYPERPOLARIZABILITY_DDD' in indices:
        if self.hyperpolarizability is None:
            ar = indices['HYPERPOLARIZABILITY_DDD']
            for ix in ar:
                tensor = zeros((3,3,3), dtype=dtype)
                indx = 0
                for ia in range(3):
                    for ib in range(3):
                        for ic in range(3):
                            tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                            if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                            indx += 1
                try:
                    self.hyperpolarizability = append(self.hyperpolarizability, tensor, axis=0)
                except ValueError:
                    self.hyperpolarizability = array([tensor])

    if 'HYPERPOLARIZABILITY_DDQ' in indices:
        ar = indices['HYPERPOLARIZABILITY_DDQ']
        for ix in ar:
            tensor = zeros((3,3,6), dtype=dtype)
            indx = 0
            for ia in range(3):
                for ib in range(3):
                    for ic in range(6):
                        tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                        if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                        indx += 1
            try:
                self.beta_ddq = append(self.beta_ddq, tensor, axis=0)
            except ValueError:
                self.beta_ddq = array([tensor])

    if 'HYPERPOLARIZABILITY_DQD' in indices:
        ar = indices['HYPERPOLARIZABILITY_DQD']
        for ix in ar:
            tensor = zeros((3,6,3), dtype=dtype)
            indx = 0
            for ia in range(3):
                for ib in range(6):
                    for ic in range(3):
                        tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                        if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                        indx += 1
            try:
                self.beta_dqd = append(self.beta_dqd, tensor, axis=0)
            except ValueError:
                self.beta_dqd = array([tensor])

    if 'HYPERPOLARIZABILITY_QDD' in indices:
        ar = indices['HYPERPOLARIZABILITY_QDD']
        for ix in ar:
            tensor = zeros((6,3,3), dtype=dtype)
            indx = 0
            for ia in range(6):
                for ib in range(3):
                    for ic in range(3):
                        tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                        if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                        indx += 1
            try:
                self.beta_qdd = append(self.beta_qdd, tensor, axis=0)
            except ValueError:
                self.beta_qdd = array([tensor])

    if 'HYPERPOLARIZABILITY_DQQ' in indices:
        ar = indices['HYPERPOLARIZABILITY_DQQ']
        for ix in ar:
            tensor = zeros((3,6,6), dtype=dtype)
            indx = 0
            for ia in range(3):
                for ib in range(6):
                    for ic in range(6):
                        tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                        if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                        indx += 1
            try:
                self.beta_dqq = append(self.beta_dqq, tensor, axis=0)
            except ValueError:
                self.beta_dqq = array([tensor])

    if 'HYPERPOLARIZABILITY_QDQ' in indices:
        ar = indices['HYPERPOLARIZABILITY_QDQ']
        for ix in ar:
            tensor = zeros((6,3,6), dtype=dtype)
            indx = 0
            for ia in range(6):
                for ib in range(3):
                    for ic in range(6):
                        tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                        if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                        indx += 1
            try:
                self.beta_qdq = append(self.beta_qdq, tensor, axis=0)
            except ValueError:
                self.beta_qdq = array([tensor])

    if 'HYPERPOLARIZABILITY_QQD' in indices:
        ar = indices['HYPERPOLARIZABILITY_QQD']
        for ix in ar:
            tensor = zeros((6,6,3), dtype=dtype)
            indx = 0
            for ia in range(6):
                for ib in range(6):
                    for ic in range(3):
                        tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                        if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                        indx += 1
            try:        
                self.beta_qqd = append(self.beta_qqd, tensor, axis=0)
            except ValueError:
                self.beta_qqd = array([tensor])

    if 'HYPERPOLARIZABILITY_QQQ' in indices:
        ar = indices['HYPERPOLARIZABILITY_QQQ']
        for ix in ar:
            tensor = zeros((6,6,6), dtype=dtype)
            indx = 0
            for ia in range(6):
                for ib in range(6):
                    for ic in range(6):
                        # Zhongwei: in case there are two many digits in a number
                        try:
                           tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[3])
                        except IndexError:
                           if "-" in f[ix+indx+ip].split()[2]:
                              tensor[ia][ib][ic] = float(f[ix+indx+ip].split()[2][2:-1])
                           else:
                              self._raise_or_pass('Error reading the beta_qqq value.')
                        if lcmplx: tensor[ia][ib][ic] += 1j*float(f[ix+indx+ip].split()[4])
                        indx += 1
            try:
                self.beta_qqq = append(self.beta_qqq, tensor, axis=0)
            except ValueError:
                self.beta_qqq = array([tensor])
