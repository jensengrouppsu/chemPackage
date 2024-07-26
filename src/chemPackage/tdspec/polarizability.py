from __future__ import print_function, division
from numpy import array, where, append, arange, zeros
from numpy import row_stack, column_stack, argsort
from chem import collect

def collect_frequency(self, f, indices):
    '''Drive the collection of polarizabilities of different types.'''

    ar = indices['NORMAL MODE']
    for ix in ar:
        mode = float(f[ix].split()[2])
        if self.v_frequencies == None:
            self.v_frequencies = array(mode)
        else:
            self.v_frequencies = append(self.v_frequencies, mode)
    self.nmodes = len(self.v_frequencies)


def collect_alpha(self, f, indices):
    '''Collects the dipole-dipole polarizability.'''

    ar = indices['POLARIZABILITY']
    for ix in ar:
        temp = zeros((1,3,3), dtype=complex)
        il = -1
        for ia in range(3):
            for ib in range(3):
                il += 1 
                temp[0][ia][ib] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.qm_pol == None:
            self.qm_pol = temp
        else:
            self.qm_pol = append(self.qm_pol, temp, axis=0)

def collect_Atensor(self, f, indices):
    '''Collects the dipole-quadrupole polarizability.'''

    ar = indices['A-TENSOR']
    for ix in ar:
        temp = zeros((1,3,6), dtype=complex)
        il = -1
        for ia in range(3):
            for ib in range(6):
                il += 1
                temp[0][ia][ib] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.atensor == None:
            self.atensor = temp
        else:
            self.atensor = append(self.atensor, temp, axis=0)

def collect_Astensor(self, f, indices):
    '''Collects the quadrupole-dipole polarizability.'''

    ar = indices['As-TENSOR']
    for ix in ar:
        temp = zeros((1,3,6), dtype=complex)
        il = -1
        for ia in range(3):
            for ib in range(6):
                il += 1
                temp[0][ia][ib] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.astensor == None:
            self.astensor = temp
        else:
            self.astensor = append(self.astensor, temp, axis=0)

def collect_Gtensor(self, f, indices):
    '''Collects the dipole-magnetic dipole polarizability.'''

    ar = indices['G-TENSOR']
    for ix in ar:
        temp = zeros((1,3,3), dtype=complex)
        il = -1
        for ia in range(3):
            for ib in range(3):
                il += 1
                temp[0][ia][ib] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.gtensor == None:
            self.gtensor = temp
        else:
            self.gtensor = append(self.gtensor, temp, axis=0)

def collect_Gstensor(self, f, indices):
    '''Collects the dipole-magnetic dipole polarizability.'''

    ar = indices['Gs-TENSOR']
    for ix in ar:
        temp = zeros((1,3,3), dtype=complex)
        il = -1
        for ia in range(3):
            for ib in range(3):
                il += 1
                temp[0][ia][ib] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.gstensor == None:
            self.gstensor = temp
        else:
            self.gstensor = append(self.gstensor, temp, axis=0)

def collect_Ctensor(self, f, indices):
    '''Collects the quadrupole-quadrupole polarizability.'''

    ar = indices['C-TENSOR']
    for ix in ar:
        temp = zeros((1,6,6), dtype=complex)
        il = -1
        for ia in range(6):
            for ib in range(6):
                il += 1
                temp[0][ia][ib] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.ctensor == None:
            self.ctensor = temp
        else:
            self.ctensor = append(self.ctensor, temp, axis=0)

def collect_Dtensor(self, f, indices):
    '''Collects the quadrupole-magnetic dipole polarizability.'''

    ar = indices['D-TENSOR']
    for ix in ar:
        temp = zeros((1,6,3), dtype=complex)
        il = -1
        for ia in range(6):
            for ib in range(3):
                il += 1
                temp[0][ia][ib] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.dtensor == None:
            self.dtensor = temp
        else:
            self.dtensor = append(self.dtensor, temp, axis=0)

def collect_Dstensor(self, f, indices):
    '''Collects the quadrupole-magnetic dipole polarizability.'''

    ar = indices['Ds-TENSOR']
    for ix in ar:
        temp = zeros((1,3,6), dtype=complex)
        il = -1
        for ia in range(6):
            for ib in range(3):
                il += 1
                temp[0][ib][ia] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.dstensor == None:
            self.dstensor = temp
        else:
            self.dstensor = append(self.dstensor, temp, axis=0)

# Zhongwei: hyperpolarizability
def collect_beta(self, f, indices):
    '''Collects the dipole-dipole-diple hyperpolarizability.'''

    ar = indices['HYPERPOLARIZABILITY']
    for ix in ar:
        temp = zeros((1,3,3,3), dtype=complex)
        il = -1
        for ia in range(3):
            for ib in range(3):
                for ic in range(3):
                    il += 1 
                    temp[0][ia][ib][ic] = float(f[ix+il].split()[1]) + float(f[ix+il].split()[2])*1.j
        if self.hyperpolarizability == None:
            self.hyperpolarizability = temp
        else:
            self.hyperpolarizability = append(self.hyperpolarizability, temp, axis=0)
