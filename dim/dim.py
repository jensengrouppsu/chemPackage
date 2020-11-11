from __future__ import print_function, division
from ..chemdata import ChemData
from ..constants import elem
from numpy import array, where, empty_like, reshape
from sys import exit

class DIM(ChemData):
    
    def __init__(self, name):
        import os

        ChemData.__init__(self)
        self.program = 'DIM'
        # Find extention
        ftype = os.path.splitext(name)[1]
        if ftype not in ('.out', '.run', '.inp', '.dim'):
            raise ValueError (ftype+' not a recognized extention')
        self.filetype = ftype[1:]
        self.filename = name

        # List of input keys.
        self.blockkeys = ('XYZ',)
        # Each element is a possible key
        for e in elem[1:]:
            self.blockkeys += (e.upper(),)
        self.linekeys = ('TOLERANCE', 'TOTALCHARGE', 'DAMPING', 'PRINTLEVEL',
                         'ALGORITHM', 'NOPRINT', 'PRINT', 'FREQRANGE',
                         'FREQUENCY', 'OUTPUT',)
        self.singlekeys = ('CPIM', 'PIM', 'NOPOL', 'NOCHAR', 'NONINTERACTING',
                           'DEBUG', 'BOHR', 'COORDDEPEND', 'DDA','COARSEGRAIN', 'RETARDATION')
    
    def _collect(self, abort=False):

        # Set the abort flag
        self._abort = abort

        # Read in file
        from .read_file import read_file
        f, indices = read_file(self)

        # Read input block
        from .input_block import collect_input
        collect_input(self, f, indices)

        # Determine calculation type
        self.__det_calc_type()

        # Return now if an input file
        if self.filetype != 'out': return

        # Techical properties
        from .dimproperties import collect_technical
        collect_technical(self, f, indices)

        # Collect all exciting information
        from .dimproperties import collect_dim
        collect_dim(self, f, indices)

        # Collect timing
        from .dimproperties import collect_timing
        collect_timing(self, f, indices)


    def __det_calc_type(self):
        '''Determine the calculation type based on what is found in the keys'''
        self.calctype.add('DIM')
        self.calctype.add('POLARIZABILITY')
        if 'PRINT' in self.key:
            if 'ATMDIP' in self.key['PRINT']:
                self.subkey.add('ATMDIP')
            if 'ENERGY' in self.key['PRINT']:
                self.subkey.add('ENERGY')
        if 'NOPRINT' in self.key:
            if 'POL' in self.key['NOPRINT']:
                self.calctype.discard('POLARIZABILITY')
            if 'ATMDIP' in self.key['NOPRINT']:
                self.subkey.discard('ATMDIP')
            if 'ENERGY' in self.key['NOPRINT']:
                self.subkey.discard('ENERGY')
        if 'CPIM' in self.key:
            self.subkey.add('CPIM')
        elif 'PIM' in self.key:
            self.subkey.add('PIM')
        if 'FREQRANGE' in self.key or 'FREQUENCY' in self.key:
            self.calctype.add('FD')
        else:
            self.calctype.add('STATIC')
