from __future__ import print_function

def part_1(psum, p, ndim):
    for n in range(ndim):
        sum = 0
        for m in range(ndim+1):
            sum += p[m][n]
        psum[n] = sum

def part_2(y):
    ilo = 0
    if (y[0] > y[1]):
        ihi = 0
        inhi = 1
    else:
        ihi = 1
        inhi = 0
    return ihi, inhi

def amotry(p, y, psum, ndim, funk, ihi, fac, opt=None):
    from numpy import zeros

    ptry = zeros((ndim), dtype=float)
    fac1 = (1.-fac)/ndim
    fac2 = fac1 - fac

    for j in range(ndim):
        ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2

    if opt is None:
        ytry = funk(ptry)
    else:
        ytry = funk(ptry, opt)

    if ytry < y[ihi]:
        y[ihi] = ytry
        for j in range(ndim):
            psum[j] = psum[j] - p[ihi][j] + ptry[j]
            p[ihi][j] = ptry[j]
    return ytry

def amoeba(p, y, ndim, funk, opt=None, ftol=1e-6, ITMAX=5000):
    '''
    A downhill simplex multi-dimension minimization
    using the Nelder and Mead method.

    Input:
    funk(x) - function to minimize where x(ndim)
              are the input variables.
    ndim - number of dimensions (variables) that funk
           accepts.
    p(ndim+1,ndim) - the ndim+1 vertices with values
                     of the ndim variables at each
                     vertex point.
    y(ndim+1) - values of funk evaluated at each of
                the ndim+1 verticies
    ftol - fractional convergence tolerance.
    ITMAX - maximum number of iterations.

    Returns:
    p - coordinates of the vertices of a minimized
        simplex.
    y - funk evaluated at each of the vertices of the
        minimized simplex.

    See pg. 404 of Numerical Recipes in FORTRAN, 1992.
    '''

    from numpy import array, zeros

    minimized = False
    call_part1 = True
    psum = zeros((ndim), dtype=float)

    niter = 0

    iter = 0
    while not minimized:
        if call_part1:
            part_1(psum, p, ndim)
            call_part1 = False
        ilo = 0
        ihi, inhi = part_2(y)
        for i in range(ndim+1):
            if y[i] < y[ilo]: ilo=i
            if y[i] > y[ihi]:
                inhi = ihi
                ihi = i
            elif y[i] > y[inhi]:
                if i != ihi: inhi = i
        rtol = 2. * abs(y[ihi]-y[ilo]) / abs(y[ihi]+y[ilo])
        niter += 1
        print ('Iteration {0:4}: tolerance = {1:10.3e}'.format(niter, rtol))

        if rtol < ftol:
            minimized = True
            swap = y[0]
            y[0] = y[ilo]
            y[ilo] = swap
            for n in range(ndim):
                swap = p[0][n]
                p[0][n] = p[ilo][n]
                p[ilo][n] = swap
            return p, y

        if (iter > ITMAX): return p, y

        iter += 2

        ytry = amotry(p, y, psum, ndim, funk, ihi, -1., opt=opt)

        if ytry < y[ilo]:
            ytry = amotry(p, y, psum, ndim, funk, ihi, 2., opt=opt)
        elif ytry >= y[inhi]:
            ysave = y[ihi]
            ytry = amotry(p, y, psum, ndim, funk, ihi, 2., opt=opt)
            if ytry > ysave:
                for i in range(ndim+1):
                    if i != ilo:
                        for j in range(ndim):
                            psum[j] = 0.5 * ( p[i][j] + p[ilo][j] )
                            p[i][j] = psum[j]
                        if opt is None:
                            y[i] = funk(psum)
                        else:
                            y[i] = funk(psum, opt)
                iter += ndim
                call_part1 = True
                continue
        else:
            iter += 1
        continue
