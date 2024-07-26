from __future__ import print_function, division
from math import hypot, atan2, sin, cos

'''Here is an assortment of functions that convert from one form to another.'''

def dhms2sec(dhms):
    '''Return the number of seconds from a number in D:H:M:S format.'''
    # Add trailing zeros if dhms ends in ':'
    if dhms[-1] == ':': dhms += '00'
    seconds = None
    nums = dhms.split(':')
    if [n for n in nums if n.isdigit()]:
        seconds = int(nums.pop()) # Last index is seconds.
        if nums: seconds += int(nums.pop()) * 60 # Next is minutes.
        if nums: seconds += int(nums.pop()) * 3600 # Next is hours.
        if nums: seconds += int(nums.pop()) * 86400 # Last is days.
    return seconds

def intround(val):
    '''Return an int from an arbitrary type, or None if not possible.'''
    try: # Easy conversion (floats)
        return int(round(val))
    except TypeError: # Complex or string.
        try:
            return int(round(float(val)))
        except TypeError: # Complex number
            try:
                return int(round(abs(val)))
            except TypeError: # Complex string
                try:
                    return int(round(abs(complex(val))))
                except TypeError:
                    return None

def string2int(string):
    '''Convert to integer if an integer string.'''
    try:
        if string.isdigit():
            return int(string)
        else:
            return string
    except:
        return string

def cart2sphere(p):
    '''Converts a point in Cartesian coordinates to spherical coordinates.
    Accepts either a three-element list or numpy array.  Returns a three
    element tuple, (r, theta, phi), where theta and phi are in radians.
    '''
    from functions import norm
    x, y, z = p[0], p[1], p[2]
    return norm(p), atan2(hypot(x, y), z), atan2(y, x)

def cart2cylinder(p):
    '''Converts a point in Cartesian coordinates to cylindrical coordinates.
    Accepts either a three-element list or numpy array.  Returns a three
    element tuple,(s, phi, z), where phi is in radians.
    '''
    x, y, z = p[0], p[1], p[2]
    return hypot(x, y), atan2(y, x), z

def sphere2cart(r, theta, phi):
    '''Converts a point in spherical coordinates to Cartesian coordinates.
    Returns a 1x3 numpy array of the coordinates. Angles must be given in
    radians.
    '''
    from numpy import array
    return array([r * sin(theta) * cos(phi),
                  r * sin(theta) * sin(phi),
                  r * cos(theta)])

def sphere2cylinder(r, theta, phi):
    '''Converts a point in spherical coordinates to cylindical coordinates.
    Returns a three element tuple, (s, phi, z). Angles must be given in
    radians and are returned in radians.
    '''
    return (hypot(r * sin(theta) * cos(phi),  r * sin(theta) * sin(phi)),
            phi,
            r * cos(theta))

def cylinder2cart(s, phi, z):
    '''Converts a point in cylindircal coordinates to Cartesian coordinates.
    Returns a 1x3 numpy array of the coordinates. phi must be given in
    radians.
    '''
    from numpy import array
    return array([s * cos(phi), s * sin(phi), z])

def cylinder2sphere(s, phi, z):
    '''Converts a point in cylindrical coordinates to spherical coordinates.
    Returns a three element tuple, (r, theta, phi). Angles must be given in
    radians and are returned in radians.
    '''
    return hypot(s, z), atan2(s, z), phi
