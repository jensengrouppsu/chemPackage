from .natsort import natsort_key, natsorted, index_natsorted
from .input_reader import InputReader, ReaderError, SUPPRESS, Namespace
from .files import *
from .compile import *

__all__ = [
           'InputReader',
           'ReaderError',
           'SUPPRESS',
           'Namespace', 
           'natsort_key',
           'natsorted',
           'index_natsorted',
           'abs_file_path',
           'locate',
           'file_safety_check',
           'range_check',
           'compile_target',
           'compile_on_fly',
          ]

def range_check(low, high, expand=False, asint=False):
    '''\
    :py:func:`range_check` will verify that that given range has a
    *low* lower than the *high*.  If both numbers are integers, it
    will return a list of the expanded range unless *expand* is
    :py:const:`False`, in which it will just return the high and low.
    If *low* or *high* is not an integers, it will return the *low*
    and *high* values as floats.

    :argument low:
       The low value if the range to check.
    :type low: float, int
    :argument high:
       The high value if the range to check.
    :type high: float, int
    :keyword expand:
        If :py:const:`True` and both *low* or *high* are integers, then
        :py:func:`range_check` will return the range of integers between
        *low* and *high*, inclusive. Otherwise, :py:func:`range_check`
        just returns *low* and *high*.
    :type expand: bool, optional
    :keyword asint:
        If *expand* is :py:const:`False`, this will attempt to return the
        *low* and *high* as integers instead of floats.
    :type expand: bool, optional
    :rtype:
        See the explanation of *expand*.
    :exception:
        * :py:exc:`ValueError`: *low* > *high*.
        * :py:exc:`ValueError`: *low* or *high* cannot be converted to a
          :py:func:`float`.
    '''
    # Convert to float first.  A ValueError is raised if not possible
    low = float(low)
    high = float(high)
    # Raise a value error if the range is invalid
    if low >= high:
        raise ValueError('low >= high')
    # If we need special integer handling, check that we have integers
    if (expand or asint) and int(low) == low and int(high) == high:
        if expand:
            return range(int(low), int(high)+1)
        else:
            return int(low), int(high)
    # Otherwise return the floats
    else:
        return low, high
