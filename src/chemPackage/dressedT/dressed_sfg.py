#! /usr/bin/env python
from __future__ import print_function

def dressed_sfg(f, E0, E1, E2, FG0, FG1, FG2, E=True, FG=True,**kwargs):

    '''Dress the beta tensor with E0(scat),E1(vis),E2(IR) and FG0, FG1, FG2 for SFG.'''
    import os
    from ..constants import BOHR2ANGSTROM as B2A, KRONECKER3 as delta, LEVICIVITA3 as epsilon, WAVENUM2HART, MASS_ELECT, AMU, DEBYE2AU
    from numpy import array, zeros, radians, cos, sin, einsum
    from dressed_func import plot, generate_field, calc_chi2_sfg, print_datafile

    #### NB: You may need to import other functions from dressed_func as needed ####
    f.collect_tensor_derivatives()
    f.dgdip = DEBYE2AU(f.dgdip)

    if E == False:
        BD = einsum('ijk,ia->ijka',f.qm_pol,f.dgdip)

    if E:
        # verify type of calculation and source of E and FG
        E0_sum = zeros(E0.shape,dtype=complex)
        E1_vis = zeros(E1.shape,dtype=complex)
        E2_ir = zeros(E2.shape,dtype=complex)

        for i in range(f.nmodes):
            E0_sum[i] = E0[i] + delta

        for i in range(f.nmodes):
            E1_vis[i] = E1[i] + delta

        for i in range(f.nmodes):
            E2_ir[i] = E2[i] + delta
       
        DDD = einsum('ijk,ia->ijka',f.qm_pol,f.dgdip)
        BD = einsum('idef,iad,ibe,icf->iabc', DDD, E0_sum, E1_vis, E2_ir)
    
    if FG:
        # calculate quadrupole-dipole dipole beta, dipole-dipole quadrupole beta, quadrupole-quadrupole dipole beta
        # quadrupole-dipole quadrupole beta and quadrupole-quadrupole quadrupole beta
        QDD = einsum('ijkl,ia->ijkla',f.atensor,f.dgdip)
        DQD = QDD
        DDQ = einsum('ijk,iab->ijkab',f.f.qm_pol,f.quadrupole)
        QQD = einsum('ijklm,ia->ijklma',f.ctensor,f.dgdip)
        QDQ = einsum('ijkl,iab->ijklab',f.atensor,f.quadrupole)
        DQQ = QDQ
        QQQ = einsum('ijklm,iab->ijklmab',f.ctensor,f.quadrupole)
        
        #print("QDD {0}".format(QDD.max()))
        #print("DDQ {0}".format(DDQ.max()))
        #print("QQD {0}".format(QQD.max()))
        #print("QDQ {0}".format(QDQ.max()))
        #print("QDQ {0}".format(QDQ.max()))
        #print("QQQ {0}".format(QQQ.max()))
        #BD += (1./3.) * einsum('idefg,iade,ibf,icg->iabc', QDD, FG0, E1_vis, E2_ir)
        #BD += (1./3.) * einsum('idefg,iad,ibef,icg->iabc', DQD, E0_sum, FG1, E2_ir)
        #BD += (1./3.) * einsum('idefg,iad,ibe,icfg->iabc', DDQ, E0_sum, E1_vis, FG2)
        #BD += (1./9.) * einsum('idefgh,iade,ibfg,ich->iabc', QQD, FG0, FG1, E2_ir)
        #BD += (1./9.) * einsum('idefgh,iade,ibf,icgh->iabc', QDQ, FG0, E1_vis, FG2)
        #BD += (1./9.) * einsum('idefgh,iad,ibef,icgh->iabc', DQQ, E0_sum, FG1, FG2)
        #BD += (1./27.) * einsum('idefghj,iade,ibfg,ichj->iabc', QQQ, FG0, FG1, FG2)
        BD += (1./3.) * einsum('iade,idefg,ibf,icg->iabc', FG0, QDD, E1_vis, E2_ir)
        BD += (1./3.) * einsum('iad,idefg,ibef,icg->iabc', E0_sum,DQD, FG1, E2_ir)
        BD += (1./3.) * einsum('iad,idefg,ibe,icfg->iabc', E0_sum,DDQ, E1_vis, FG2)
        BD += (1./9.) * einsum('iade,idefgh,ibfg,ich->iabc', FG0, QQD, FG1, E2_ir)
        BD += (1./9.) * einsum('iade,idefgh,ibf,icgh->iabc', FG0, QDQ, E1_vis, FG2)
        BD += (1./9.) * einsum('iad,idefgh,ibef,icgh->iabc', E0_sum, DQQ, FG1, FG2)
        BD += (1./27.) * einsum('iade,idefghj,ibfg,ichj->iabc', FG0, QQQ, FG1, FG2)

    pre = -1/(2*WAVENUM2HART(f.v_frequencies))
    omega_ir = WAVENUM2HART(f.v_frequencies)
    gam = 10j
    temp = 1 / (WAVENUM2HART(f.v_frequencies)-WAVENUM2HART(f.v_frequencies)+WAVENUM2HART(gam))                             
    beta = einsum('i,iabc,i->iabc',pre,BD,temp)
    return beta

def return_kwargs(string, val=False, **kwargs):
    '''Check if a string is contained within the given kwargs
    and return the value of that kwargs, otherwise return val.'''
    if string in kwargs:
        return kwargs[string]
    else:
        return val

