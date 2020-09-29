from __future__ import print_function, division
from chemdata import ChemData
from errorclass import CollectionError
from numpy import array, append
import itertools

class CubeFile(ChemData):
    '''Read and contain data from Gaussian cube files.'''

    # Instantiate the CubeFile class.
    def __init__(self, filename):
        '''Set up the CubeFile class.'''
        import os

        ChemData.__init__(self)
        ftype = os.path.splitext(filename)[1]
        self.filetype = ftype[1:]
        self.filename = filename

    def _collect(self, abort=False):
        '''Collect data from a Gaussian cube file.'''  
         
        # Save the abort status
        self._abort = abort

        # Collect all data into memory for data retention
        with open(self.filename) as fl:
            f = tuple([line.rstrip() for line in fl]) 

        # Collect the total number of atoms.        
        self.natoms = int(f[2].split()[0])

        # Collect the title of the cube file (if it exists).
        if f[0] == '':
            self.title = '\n'
        else:
            self.title = f[0]

        # Collect the "dimension" of the MO.  This allows the storage
        # of the MO with the same dimensions as the file made by the 
        # quantum chemistry software.
        self.orbitaldim = array([], dtype=int)
        for ln in f[3:6]:
            ln = ln.split()[0]
            self.orbitaldim = append(self.orbitaldim, int(ln))

        # Because of how we collect the MO, it is first a list
        # with self.orbitaldim[0]*self.orbitaldim[1]*self.orbitaldim[2]
        # elements.  To stay consistent with the original cube file,
        # we need to change the shape of the list.  The numpy array 
        # containing the MO itself will have the dimension:
        #
        # self.orbitaldim[1]*self.orbitaldim[2],self.orbitaldim[0]
 
        start = self.natoms + 6
        MO = itertools.chain.from_iterable([x.split() for i,x in enumerate(f[start:])])
        MO = array(list(MO),dtype=float)
        MO.shape = (self.orbitaldim[1]*self.orbitaldim[2],self.orbitaldim[0])
        self.orbitalplot = MO

