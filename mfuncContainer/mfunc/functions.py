from __future__ import print_function, division
from math import hypot, atan2, sin, cos, sqrt, fmod, pi
from mfunc.mvmult import vecvecr, vecvecc, matmatr, matmatc
from mfunc.mvmult import vecmatr, vecmatc, matvecr, matvecc
from mfunc.mvmult import normr, normc

'''Here is an assortment of math functions.'''

__all__ = [
           'lorentzian',
           'sum_lorentzian',
           'gaussian',
           'sum_gaussian',
           'vectorangle',
           'norm',
           'rms',
           'quadrant',
           'mult',
           'roundint',
          ]

def mult(x, y):
    '''Performs vector-vector multiplication (dot product),
    matrix-vector multiplication, or matrix-matrix multiplication.
    It is a wrapper for MKL routines.  Automatically determines if
    you need level 1, 2 or 3 BLAS. 

    IF YOUR ARRAYS ARE IN C-ORDER:
    mult will always be faster than np.dot (dot) for vector-vector
    multiplication.  However, it will only be faster than dot for
    large matrix-vector multiplications or moderate matrix-matrix
    multiplications.  The reason is that for matrices (not vectors)
    the array must be reordered in memory before being passed to the
    FORTRAN code.  The slows down the initiallization matrix-vector
    and matrix-matrix routines. Matrix-vector multiplication incurs
    a more severe penalty (counterintuitively, since there is only one
    matrix, not two) because dot and MKL routine are comparable in
    speed, whereas MKL's matrix-matrix multiplication is much faster
    than dot.

    IF YOUR ARRAYS ARE IN FORTRAN-ORDER:
    mult will always be faster than dot (more than twice as fast
    for as small as 50x50 * 50x50 multiplication).  This is because
    there is no initiallization penalty since no reordering must be
    done.
    '''
    from numpy import asarray
    x = asarray(x)
    y = asarray(y)
    # Vector-vector multiplication
    if x.ndim == 1 and y.ndim == 1:
        # Make sure that the vector lengths are the same
        if x.size != y.size:
            raise ValueError('mult(): Vector length mismatch, '
                             '({0} != {1})'.format(x.size, y.size))
        else:
            if 'complex' in str(x.dtype) or 'complex' in str(y.dtype):
                return vecvecc(x, y)
            else:
                return vecvecr(x, y)
    # Matrix-vector multiplication
    elif x.ndim > 1 and y.ndim ==1:
        # Verify that vector lenths are correct
        if x.shape[1] != y.size:
            raise ValueError('mult(): Vector length must match matrix col '
                             '({0} != {1})'.format(y.size, x.shape[1]))
        else:
            if 'complex' in str(x.dtype) or 'complex' in str(y.dtype):
                return matvecc(x, y)
            else:
                return matvecr(x, y)
    # Vector-matrix multiplication
    elif x.ndim == 1 and y.ndim > 1:
        # Verify that vector lenths are correct
        if y.shape[0] != x.size:
            raise ValueError('mult(): Matrix row must match vector length '
                             '({0} != {1})'.format(x.size, y.shape[0]))
        else:
            if 'complex' in str(x.dtype) or 'complex' in str(y.dtype):
                return vecmatc(x, y)
            else:
                return vecmatr(x, y)
    # Matrix-matrix multiplication
    else:
        # Verify that the inner dimentions match
        if x.shape[1] != y.shape[0]:
            raise ValueError('mult(): A col must match B row '
                             '({0} != {1})'.format(x.shape[1], y.shape[0]))
        else:
            if 'complex' in str(x.dtype) or 'complex' in str(y.dtype):
                return matmatc(x, y)
            else:
                return matmatr(x, y)

def norm(x):
    '''Takes the norm of x, i.e. sqrt(dot(x, x)).  Accepts float or
    complex, but only returns float.  It is a wrapper around MKL's
    dnrm2 and dznrm2.  Can accept a list or numpy array.
    '''
    from numpy import asarray
    x = asarray(x)
    if 'complex' in str(x.dtype):
        return normc(x)
    else:
        return normr(x)

def lorentzian(x, peak=0, height=1.0, fwhm=None, hwhm=None):
    '''Calculates a three-parameter lorentzian for a given domain.'''
    if fwhm is not None and hwhm is not None:
        raise ValueError ('lorentzian: Onle one of fwhm or hwhm must be given')
    elif fwhm is not None:
        gamma = fwhm / 2
    elif hwhm is not None:
        gamma = hwhm
    else:
        gamma = 0.1
    # pi is included as a normalization factor
    return  ( height / pi ) * ( gamma / ( ( x - peak )**2 + gamma**2 ) )

def sum_lorentzian(x, peak=None, height=None, fwhm=None, hwhm=None):
    '''Calculates and sums several lorentzians to make a spectrum.

    'peak' and 'height' are numpy arrays of the peaks and heights that
    each component lorentzian has.
    '''
    from numpy import array
    if peak is None or height is None:
        raise ValueError ('Must pass in values for peak and height')
    if peak.shape != height.shape:
        raise ValueError ('peak and height must be the same shape')

    l = lorentzian
    y = array([l(x,peak[i], height[i], fwhm, hwhm) for i in range(len(peak))])
    return y.sum(axis=0)

def gaussian(x, peak=0, height=1.0, fwhm=None, hwhm=None):
    '''Calculates a three-parameter gaussian for a given domain.'''
    import numpy as np
    if fwhm is not None and hwhm is not None:
        raise ValueError ('gaussian: Onle one of fwhm or hwhm must be given')
    elif fwhm is not None:
        gamma = fwhm / 2
    elif hwhm is not None:
        gamma = hwhm
    else:
        gamma = 0.1
    # use np.log and np.exp to avoid errors
    sigma = gamma / sqrt( np.log(4) )
    # pi is included as a normalization factor
    return  ( height / ( sigma * sqrt(2 * pi) ) ) * np.exp( -(x - peak)**2 / ( 2 * sigma**2 ) )

def sum_gaussian(x, peak=None, height=None, fwhm=None, hwhm=None):
    '''Calculates and sums several gaussians to make a spectrum.

    'peak' and 'height' are numpy arrays of the peaks and heights that
    each component gaussian has.
    '''
    from numpy import array
    if peak is None or height is None:
        raise ValueError ('Must pass in values for peak and height')
    if peak.shape != height.shape:
        raise ValueError ('peak and height must be the same shape')

    l = gaussian
    y = array([l(x,peak[i], height[i], fwhm, hwhm) for i in range(len(peak))])
    return y.sum(axis=0)

def vectorangle(v1, v2):
    '''Determine the angle between two vectors in radians.
    Both vectors must be of same dimension.
    '''
    from numpy import cross
    if len(v1) != len(v2):
        raise ValueError ('vectorangle(): len(v1) != len(v2)')
    elif mult(v1, v1) == 0 or mult(v2, v2) == 0:
        raise ValueError ('vectorangle(): Zero-length vector')

    # phi = atan(||v1 x v2||, v1 . v2)
    return atan2(norm(cross(v1, v2)), mult(v1, v2))

def rms(vec):
    '''Calculates the root mean square of a vector.'''
    return sqrt(mult(vec, vec) / len(vec))

def roundint(num, integer=1):
    '''Rounds the given number to the nearest multiple of the given integer.
    If the given number is a float, it is first cast as an integer.'''
    integer = int(integer)
    return int(integer * round(num / integer))

def quadrant(angle):
    '''Determines what quadrant an angle is in.  If the angle is on an axis,
    it returns the quadrant above it, i.e. 90 degrees returns 2, not 1.
    Angle must be given in radians.
    '''
    reduced = fmod(angle, 360)
    s = round(sin(reduced),14)
    c = round(cos(reduced),14)
    if s >= 0 and c > 0:
        return 1
    elif s > 0 and c <= 0:
        return 2
    elif s <=0 and c < 0:
        return 3
    elif s < 0 and c >= 0:
        return 4
