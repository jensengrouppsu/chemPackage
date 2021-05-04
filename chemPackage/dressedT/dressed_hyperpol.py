#! /usr/bin/env python
from __future__ import print_function

def dressed_hyperpol(f, d=None, E=None, com=None, r=None, debug=False, e=None,
    SFG=False, hyperR=False, freq_VIS=354, **kwargs):

    from ..pol import Polarizability 
    import os
    from ..constants import BOHR2ANGSTROM as B2A, KRONECKER3 as delta, LEVICIVITA3 as epsilon
    from numpy import array, zeros, radians, cos, sin
    from dressed_func import plot, generate_field, calc_chi2_sfg
    #### NB: You may need to import other functions from dressed_func as needed ####

    if hyperR:
        infile = open('garbage','r')
        lines= infile.readlines()
        modes = []
        a=0

        while a < len(lines):
            temp = "".join(lines[a].split())
            modes.append(temp[:7])
            a = a + 10

        hpol = zeros([len(modes),3,3,3])

        c = 0
        a = 1
        while a < len(lines):
            temp = zeros([9,3])
            for n in range(9):
                temp2 = lines[a+n].split()
                for bb in range(3):
                    temp[n][bb] = temp2[bb]

            line = 0
            for i in range(3):
                for j in range(3):
                    hpol[c,i,j,0] = temp[line][0]
                    hpol[c,i,j,1] = temp[line][1]
                    hpol[c,i,j,2] = temp[line][2]
                    line += 1

            a = a + 10
            c = c +1
        beta = zeros([len(modes)/2,3,3,3], dtype=complex)
        for n in range(len(modes)/2):
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        beta[n][i][j][k] = hpol[n][i][j][k] + hpol[n+len(modes)/2][i][j][k]*1j

    # initialize x and y arrays: x - wavenumber; y - intensity
    if hyperR:
        f.v_frequencies = modes[:len(modes)/2]
        f.v_frequencies = [float(n) for n in f.v_frequencies]
    x = array(f.v_frequencies)
    y = zeros((len(x)),dtype=float)
    E0, E2 = [zeros((len(f.v_frequencies),3,3), dtype=int) for n in range(2)]
    E1 = zeros((3,3))

    # verify type of calculation and source of E and FG
    if r != None:
        if debug: print ('Calculating E and FG from an isotropic sphere...')
        E0, E1,E2 = generate_field(f, r=r, cm=com, alpha=None, e=e, SFG=SFG, hyperR= hyperR, freq_VIS = freq_VIS)
    
    if SFG:
        gam = 10j
        # cycle through all modes
        beta = zeros((len(x),3,3,3), dtype=complex)
        for n in range(len(x)):
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        temp = 1 / (-gam)
                        beta[n][i][j][k] = beta[n][i][j][k] + temp*f.polarizability[n][i][j]*f.dgdip[n][k]
    theta = 158
    psi = 0
    beta_orig =beta
#    for t in range(0,180,10):
#        for p in range(0,360,10):
#    print(t,p)
#    theta = t
#    psi = p
    beta = beta_orig
    theta = radians(theta)
    psi = radians(psi)
    phi = radians(0)
    #rotate beta
    a11 = - sin(psi)*sin(phi) + cos(theta)*cos(psi)*cos(phi)
    a12 = -cos(psi)*sin(phi) - cos(theta)*sin(psi)*cos(phi)
    a13 = sin(theta)*cos(phi)
    a21 = sin(psi)*cos(phi) + cos(theta)*cos(psi)*sin(phi)
    a22 = cos(psi)*cos(phi) - cos(theta)*sin(psi)*sin(phi)
    a23 = sin(theta)*sin(phi)
    a31 = -sin(theta)*cos(psi)
    a32 = sin(theta)*sin(psi)
    a33 = cos(theta)
    
    
    r = array([[ a11, a12, a13 ],
               [ a21, a22, a23 ],
               [ a31, a32, a33 ]], dtype=float)
    
    b = zeros((len(x),3,3,3), dtype=complex)
    for n in range(len(x)):
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for o in range(3):
                                b[n][i][j][k] =  b[n][i][j][k] + beta[n][l][m][o]*r[i][l]*r[j][m]*r[k][o]
    beta = b 
    
    b = zeros((len(x),3,3,3), dtype=complex)
    for n in range(len(x)):
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for o in range(3):
                                if SFG:
                                    b[n][i][j][k] += (E0[n][i][l] + delta[i][l]) * beta[n][l][m][o] * (E1[j][m] + delta[j][m])*(E2[n][k][o] +delta[k][o])
                                elif hyperR:
                                     b[n][i][j][k] += (E0[n][i][l] + delta[i][l]) * beta[n][l][m][o] * (E1[j][m] + delta[j][m])*(E1[k][o] +delta[k][o])
                                 
    beta = b
    if SFG:
        chi2 = calc_chi2_sfg(f, x,beta)
    
    hyperR = True
    if hyperR:
        f.calctype.add('HYPERPOLARIZABILITY')
        f.dhpol = beta
        f.nmodes = len(f.v_frequencies)
        intens = f.hpol_average()
#         outfile = ('t_'+ str(t) + 'p_' + str(p))
#         out=open(outfile, 'w')
        for i in range(f.nmodes):
            print(x[i],intens[i])
#            out.write(str(x[i])+'\t'+str(intens[i])+'\n')
#        out.close()
        #plot(Freq=x, Intensity=intens)
