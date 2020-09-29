from __future__ import print_function, division
from chemdata import ChemData
from errorclass import CollectionError
from numpy import array, append
import itertools

class EmptyFile(ChemData):
    '''Keep track of files that are empty.'''

    # Instantiate the EmptyFile class.
    def __init__(self, filename):
        '''Set up the EmptyFile class.'''
        import os

        ChemData.__init__(self)
        ftype = os.path.splitext(filename)[1]
        self.filetype = ftype[1:]
        self.filename = filename

    def _collect(self, abort=False):
        '''Store information about an Empty file.'''  
         
        # Save the abort status
        self._abort = abort

        # Calculation type is set to empty
        self.calctype.add('EMPTY')

        # Add termination status
        self.termination = 'File is empty'
