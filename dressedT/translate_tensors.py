#! /usr/bin/env python

def translate_polarizabilities(obj, R=None):
    '''Translates the A, G, D and D tensors from a ChemData obj
    using vector R (bohrs).
    if R is not given, assumes the object's center of
    nuclear charge as the vector to translate by.'''

    from ..constants import ANGSTROM2BOHR as A2B
    from copy import deepcopy
    from numpy import array

    # Make a copy of the ChemData object
    objT = deepcopy(obj)

    # Figure out how much to translate by
    if R is None:
        R = array(obj.center_of_nuc_charge * A2B, dtype=float)

    # Translate tensors
    objT.atensor = translate_Atensor(obj, R)
    objT.gtensor = translate_Gtensor(obj, R)
    objT.ctensor = translate_Ctensor(obj, R)
    objT.dtensor = translate_Dtensor(obj, R)

    # Translate coordinates
    objT.coordinates = objT.coordinates - R / A2B

    # Return translated ChemData obj
    return objT

def translate_multipoles(obj, R=None):
    '''Translates the multipole moments.'''

    from numpy import array, zeros
    from ..constants import KRONECKER3 as dt
    from ..constants import LEVICIVITA3 as eps
    from ..constants import WAVENUM2HART as F2AU
    from ..constants import ANGSTROM2BOHR as A2B
    from copy import deepcopy

    # Make a copy of the ChemData object
    objT = deepcopy(obj)

    # Figure out how much to translate by
    if R is None:
        R = obj.center_of_nuc_charge * A2B

    # Translate the quadrupole moment
    if obj.quadrupole is not None and obj.dipole is not None:
        q = deepcopy(obj.quadrupole)
        u = deepcopy(obj.dipole)
        objT.quadrupole = zeros((6), dtype=float)
        for ic in range(6):
            ia, ib = (ic*2/5), (ic/3)+(ic%3)-(ic/5)
#            ic = (ia+1)*(ib+1) - (ia+1)*(ib+1) * (ia+1)*(ib+1)/5 * (2*(ia+1)*(ib+1) - 6)/4 - 1
            objT.quadrupole[ic] += q[ic] + 1.5 * R[ib] * u[ia] + 1.5 * R[ia] * u[ib]
            if ia == ib:
                objT.quadrupole[ic] -= ( R[0]*u[0] + R[1]*u[1] + R[2]*u[2] )

    # Return translated ChemData obj
    return objT

def translate_Atensor(obj, R=[0,0,0]):
    '''Translate the A-tensor.'''

    from copy import deepcopy
    from ..constants import KRONECKER3 as dt
    from numpy import einsum, zeros

    if obj.atensor is not None and obj.polarizability is not None:
        obj.fold_Atensor('UNFOLD')
        A = deepcopy(obj.atensor)
        A -= 1.5 * einsum('c,mab->mabc', R, obj.polarizability)
        A -= 1.5 * einsum('b,mac->mabc', R, obj.polarizability)
        A += einsum('bc,d,mad->mabc', dt, R, obj.polarizability)
        return A
    else:
        return obj.atensor

def translate_Gtensor(obj, R=[0,0,0]):
    '''Translate the G-tensor.'''

    from copy import deepcopy
    from ..constants import LEVICIVITA3 as eps
    from numpy import einsum

    if obj.gtensor is not None and obj.polarizability is not None:
        G = deepcopy(obj.gtensor)
        G += 0.5 * einsum('m,bcd,c,mad->mab', obj.e_frequencies, eps, R, obj.polarizability)
        return G
    else:
        return obj.gtensor

def translate_Ctensor(obj, R=[0,0,0]):
    '''Translate the C-tensor.'''

    from copy import deepcopy
    from ..constants import KRONECKER3 as dt
    from numpy import einsum

    if obj.ctensor is not None and obj.atensor is not None and obj.polarizability is not None:
        obj.fold_Ctensor('UNFOLD')
        obj.fold_Atensor('UNFOLD')
        C = deepcopy(obj.ctensor)
        C -= 1.5 * einsum('b,macd->mabcd', R, obj.atensor)
        C -= 1.5 * einsum('a,mbcd->mabcd', R, obj.atensor)
        C -= 1.5 * einsum('d,mcab->mabcd', R, obj.atensor)
        C -= 1.5 * einsum('c,mdab->mabcd', R, obj.atensor)
        C += 2.25 * einsum('d,b,mac->mabcd', R, R, obj.polarizability)
        C += 2.25 * einsum('d,a,mbc->mabcd', R, R, obj.polarizability)
        C += 2.25 * einsum('c,b,mad->mabcd', R, R, obj.polarizability)
        C += 2.25 * einsum('c,a,mbd->mabcd', R, R, obj.polarizability)
        C += einsum('ab,e,mecd->mabcd', dt, R, obj.atensor)
        C += einsum('cd,e,meab->mabcd', dt, R, obj.atensor)
        C -= 1.5 * einsum('ab,d,e,mec->mabcd', dt, R, R, obj.polarizability)
        C -= 1.5 * einsum('ab,c,e,med->mabcd', dt, R, R, obj.polarizability)
        C -= 1.5 * einsum('cd,e,b,mae->mabcd', dt, R, R, obj.polarizability)
        C -= 1.5 * einsum('cd,e,a,mbe->mabcd', dt, R, R, obj.polarizability)
        C += einsum('ab,cd,e,g,meg->mabcd', dt, dt, R, R, obj.polarizability)
        return C
    else:
        return obj.ctensor

def translate_Dtensor(obj, R=[0,0,0]):
    '''Translate the D-tensor.'''

    from copy import deepcopy
    from ..constants import KRONECKER3 as dt
    from ..constants import LEVICIVITA3 as eps
    from numpy import einsum

    if obj.dtensor is not None and obj.polarizability is not None and obj.gtensor is not None and obj.atensor is not None:
        obj.fold_Dtensor('UNFOLD')
        obj.fold_Atensor('UNFOLD')
        D = deepcopy(obj.dtensor)
        D -= 1.5 * einsum('b,mac->mabc', R, obj.gtensor)
        D -= 1.5 * einsum('a,mbc->mabc', R, obj.gtensor)
        D += einsum('ab,d,mdc->mabc', dt, R, obj.gtensor)
        D += 0.5 * einsum('m,ced,e,mdab->mabc', obj.e_frequencies, eps, R, obj.atensor)
        D -= 0.75 * einsum('m,ced,b,e,mad->mabc', obj.e_frequencies, eps, R, R, obj.polarizability)
        D -= 0.75 * einsum('m,ced,a,e,mbd->mabc', obj.e_frequencies, eps, R, R, obj.polarizability)
        D += 0.5 * einsum('m,ab,ceg,d,e,mdg->mabc', obj.e_frequencies, dt, eps, R, R, obj.polarizability)
        return D
    else:
        return obj.dtensor
