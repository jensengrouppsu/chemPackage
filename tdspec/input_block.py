from __future__ import print_function, division
from numpy import array
import os

def collect_input(self, f, indices):
    '''Collect from the input section of the file.'''

    # Very simple collection at the moment since the TDSPEC input file
    # is not complicated.
    for item in indices.keys():
        if item == 'OPA START':
            self.key['OPA'] = True
        elif item == 'RRS START':
            self.key['RRS'] = True
        elif item == 'TPA START':
            self.key['TPA'] = True
        elif item == 'RHRS START':
            self.key['RHRS'] = True
