from __future__ import print_function, division
from math import hypot, atan2, sin, cos, sqrt, fmod, pi

'''Here is an assortment of math functions.'''

__all__ = [
           'lorentzian',
           'sum_lorentzian',
           'gaussian',
           'sum_gaussian',
           'quadrant',
           'roundint',
          ]


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
