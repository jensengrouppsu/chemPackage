from mfunc.conversions import *
from mfunc.fit import *
try:
    from mfunc.functions import *
except ImportError:
    from prep import compile_on_fly
    compile_on_fly('mfunc', 'mvmult')
    from mfunc.functions import *
try:
    from mfunc.erf import error_function as erf
except ImportError:
    from prep import compile_on_fly
    compile_on_fly('mfunc', 'erf')
    from mfunc.erf import error_function as erf

__all__ = [ 
           'fit_func',
           'fit',
           'Parameter',
           'tick_scale',
           'lorentzian',
           'sum_lorentzian',
           'gaussian',
           'sum_gaussian',
           'vectorangle',
           'mult',
           'norm',
           'rms',
           'roundint',
           'erf',
           'quadrant',
           'dhms2sec',
           'intround',
           'string2int',
           'cart2sphere',
           'cart2cylinder',
           'sphere2cart',
           'sphere2cylinder',
           'cylinder2sphere',
           'cylinder2cart',
          ]

