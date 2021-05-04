from __future__ import print_function, division
from . import constants
# Not every platform has the backends necessary for mpl, or has vpython
#try:
#    import drawing
#except ImportError:
#    pass
#import F2drawing
from .chemdata import ChemData
from .errorclass import *
from .cheminfo import __doc__
import sys
from traceback import print_exc

__all__ = [ 'collect',
            'constants',
            'ChemData',
            'ChemDataError',
            'CollectionError'
          ]

def collect(file, raise_err=True, project='all'):
    '''
    Function to determine what type of file is given and collect
    appropriately.
    
    :param file: Path to file.  Accepts relative or absolute paths.
    :type file: File Object or Pickle
    :param raise_err: Raise :py:exc:`CollectionError` on error.  Otherwise,
    execution continues.
    :rtype: ChemObj

    This is essentially a front-end for the classes that hold data from
    files.  It determines what filetype has been passed in and chooses
    the correct class.

    This will raise :py:exc:`IOError` if the file does not exist.

    On error, the default is to raise a CollectionError with complete
    traceback.  If raise_err, it will attempt to complete collection
    and ignore errors.

    '''
    import sys, os
    from .errorclass import CollectionError
    import stat

    # Check if this is a pickle.  If so, we may skip the rest of this function
    #print( os.path.splitext(file))
    if os.path.splitext(file)[1] == '.pickle':
        from cPickle import load
        # Raises IOError if doesn't exist.
        return load(open(file, 'rb'))

    # Perform tests up front.  Load file into memory to search fast
    # Raises IOError if doesn't exist.
    f = open(file).read()
    # Check if a keyword exists
    adf_out = 'Amsterdam Density Functional' in f
    adf_in = r'$ADFBIN/adf' in f
    nwchem_out = 'Northwest Computational Chemistry Package' in f
    nwchem_in = os.path.splitext(file)[1] == '.nw'
    dalton_out = 'DALTON' in f
    dalton_in = '**DALTON INPUT' in f
    qchem_out = 'Q-Chem' in f
    qchem_in = '$rem' in f or '$REM' in f
    gauss_out = 'Gaussian' in f
    gauss_in = os.path.splitext(file)[1] == '.g09'
    tdspec_out = 'TDSPEC' in f
    dim_out = 'D I S C R E T E   I N T E R A C T I O N   M O D E L' in f
    dim_in = os.path.splitext(file)[1] == '.dim'
    # In XYZ, first line is number stating # of atoms.
    xyz_file = f[0:f.find('\n')].strip().isdigit()
    # Sanely named cube files have a .cube extension.
    cube_file = os.path.splitext(file)[1] == '.cube' 
    # Check if a file is empty, which can happen sometimes due to errors
    # on the cluster.
    empty_file = os.stat(file).st_size == 0 
    del f

    # Create the ChemData object
    try:
        if nwchem_out or nwchem_in:
            from .nwchem import NWChem
            d = NWChem(file)
        elif adf_out or adf_in:
            from .adf import ADF
            d = ADF(file)
        elif dim_out or dim_in:
            from .dim import DIM
            d = DIM(file)
        elif xyz_file:
            from .chemdata import ChemData
            d = ChemData(file)
        elif tdspec_out:
            from .tdspec import TDSPEC
            d = TDSPEC(file)
        elif empty_file:
            from emptyfile import EmptyFile
            d = EmptyFile(file)            
        else:
            raise NotImplementedError
    # Invalid type extention
    except ValueError as v:
        if raise_err:
            print_exc()
            raise CollectionError('Possibly a bad extention: ', file)
        else:
            return None
    # Unknown filetype
    except NotImplementedError:
        if raise_err:
            raise CollectionError('Filetype not yet implemented: '+file)
        else:
            print('Skipping', file+'...', file=sys.stderr)
            return None

    # Collect from file
    try:
        d._collect(abort=raise_err)
    # If a known error occured  #xing
    except CollectionError as c:
        if raise_err:
            raise CollectionError (str(c) + ': ' + file)
        else:
            d.termination = c.msg
            return d
    # If a known unknown error occured.  
    except (IndexError, StopIteration, KeyError):
        msg = 'Did the calculation end prematurely?'
        if raise_err:
            print('/\\'*int(len(msg)/2), file=sys.stderr)
            print(msg, file=sys.stderr)
            print('/\\'*int(len(msg)/2), file=sys.stderr)
            print(file=sys.stderr)
            raise
        else:
            d.termination = msg
            return d
    except ValueError:
        msg = 'Possibly a number too large for the format (********)\n'
        msg += 'or a string and number running together (WORD1.785)'
        if raise_err:
            print('/\\'*27, file=sys.stderr)
            print(msg, file=sys.stderr)
            print('/\\'*27, file=sys.stderr)
            print(file=sys.stderr)
            raise
        else:
            d.termination = msg
            return d
    # Keyboard interrupt
    except KeyboardInterrupt:
        sys.exit(1)
    # Everything actually went well
    else:
        return d
