def rotation_matrix_3D(a, b, c, convention='RxRzRx', angles='deg'):
    '''Determines the rotation matrix from 3 Euler angles
    using a specified rotation convention.'''

    from numpy import array
    from math import sin, cos, pi

    # convert to radians
    if angles == 'deg':
        a = a * pi / 180.
        b = b * pi / 180.
        c = c * pi / 180.

    # find sin and cos of each angle
    s1 = sin(a)
    c1 = cos(a)
    s2 = sin(b)
    c2 = cos(b)
    s3 = sin(c)
    c3 = cos(c)

    # generate rotation matrix according to specified conventions
    if convention == 'RxRzRx':
        R = array([[ c2   , -c3*s2        , s2*s3           ],
                   [ c1*s2, c1*c2*c3-s1*s3, -c3*s1-c1*c2*s3 ],
                   [ s1*s2, c1*s3+c2*c3*s1, c1*c3-c2*s1*s3  ]])
    elif convention == 'RxRzRy':
        R = array([[ c2*c3         , -s2  , c2*s3          ],
                   [ s1*s3+c1*c3*s2, c1*c2, c1*s2*s3-c3*s1 ],
                   [ c3*s1*s2-c1*s3, c2*s1, c1*c3+s1*s2*s3 ]])
    elif convention == 'RxRyRx':
        R = array([[ c2    , s2*s3         , c3*s2           ],
                   [ s1*s2 , c1*c3-c2*s1*s3, -c1*s3-c2*c3*s1 ],
                   [ -c1*s2, c3*s1+c1*c2*s3, c1*c2*c3-s1*s3  ]])
    elif convention == 'RzRxRz':
        R = array([[ c1*c3-c2*s1*s3, -c1*s3-c2*c3*s1, s1*s2  ],
                   [ c3*s1+c1*c2*s3, c1*c2*c3-s1*s3 , -c1*s2 ],
                   [ s2*s3         , c3*s2          , c2     ]])
    # Zhongwei: two-angle rotations, copied from "dressed_hypol.py"
    #           a -> theta, b-> psi, c-> phi = 0
    elif convention == 'twoAngle':
        R = array([[ -s2*s3+c1*c2*c3, -c2*s3-c1*s2*c3, s1*c3 ],
                   [ s2*c3+c1*c2*s3 , c2*c3-c1*s2*s3 , s1*s3 ],
                   [ -s1*c2         , s1*s2          , c1    ]])
    else:
        print ('Convention unknown! Using default: RxRzRx')
        R = array([[ c2   , -c3*s2        , s2*s3           ],
                   [ c1*s2, c1*c2*c3-s1*s3, -c3*s1-c1*c2*s3 ],
                   [ s1*s2, c1*s3+c2*c3*s1, c1*c3-c2*s1*s3  ]])

    # return rotation matrix
    return R

def rotate3Dtensor(t, R=None, theta=0, phi=0, chi=0, rank=None, angles='deg', convention='RxRzRx'):
    '''An operation to rotate a tensor about Euler angles
    theta, phi and chi or rotation matrix R.
    Angles may be given in degrees (deg) or radians (rad).'''

    from numpy import einsum

    # determine the tensor rank
    if rank is None:
        rank = t.ndim

    # generate rotation matrix
    if R is None:
        R = rotation_matrix_3D(theta, phi, chi, convention=convention, angles=angles)

    # rotate the tensor
    #if rank == 1:
    #    T = einsum('ij,j', R, t)
    #elif rank == 2:
    #    T = einsum('ik,jm,km', R, R, t)
    #elif rank == 3:
    #    T = einsum('il,jm,kn', R, R, R, t)
    #elif rank == 4:
    #    T = einsum('ai,bj,ck,dl,ijkl', R, R, R, R, t)
    # Zhongwei: use the way Phil did for two-angle rotation
    if rank == 1:
        T = einsum('i,mi->m', t, R)
    elif rank == 2:
        T = einsum('ij,mj->im', t, R)
    elif rank == 3:
        T = einsum('ijk,mj,nk->imn', t, R, R)
    elif rank == 4:
        T = einsum('ijkl,mj,nk,ol->imno', t, R, R, R)

    return T
