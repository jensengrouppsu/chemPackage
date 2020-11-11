from __future__ import print_function
from numpy import array, zeros, arange

class density_class():
    '''Class to calculate, collect and read the electron density.
    The electron density is calculated in :py:attr:`TAPE41`,
    by using the :py:attr:`$ADFBIN/cjdensf` program.'''

    def __init__(self, filename, type='SCF', grid='Medium', extend=2.0,
                 t21=None, t41=None, frags=False):

        import os
        import pyadf
        from subprocess import call

        # Search for TAPE41 or TAPE21 in the samedirectory
        cwd = os.getcwd()
        if not t21:
            t21 = filename[:-4]+'.t21'
        if not t41:
            t41 = filename[:-4]+'.t41'

        if not (os.path.isfile(t41)):
            print ('TAPE41 not found, calculating from TAPE21.')
            if (os.path.isfile(t21)):

                # Calculate TAPE41 from TAPE21
                str = 'cp {0} TAPE21'.format(t21)
                er = call(str, shell=True)
                fragstr     = 'Fragment Active'
                if (frags): fragstr = ' AllFragments'
                denstr      = 'Dens {0}'.format(type)
                extendstr   = 'Extend {0}'.format(extend)
                gridstr     = 'Grid Save {0}'.format(grid)
                str = '$ADFBIN/cjdensf << eor > {0}.cjdensf\n'  \
                      ' {1}\n'                         \
                      ' {2}\n'                                  \
                      ' {3}\n'                                  \
                      ' {4}\n'                                  \
                      'eor'.format(filename[:-4], fragstr, denstr, extendstr, gridstr)
                er = call(str, shell=True)
                er = call('rm TAPE21', shell=True)
                str = 'mv TAPE41 {0}'.format(t41)
                er = call(str, shell=True)
                str = 'cat logfile >> {0}'.format(filename[:-4]+'.cjdensf')
                er = call(str, shell=True)
                er = call('rm logfile', shell=True)

            else:
                print ('TAPE21 not found. Density Cannot be calculated!')
                return

        # Read geometry from TAPE41
        obj = pyadf.kf.kffile(t41)
        self.nnuc = obj.read('Geometry', 'nnuc')[0]
        self.qtch = obj.read('Geometry', 'qtch')

        # Get grid coordinates
        self.start_point    = obj.read('Grid', 'Start_point')
        self.nx             = obj.read('Grid', 'nr of points x')[0]
        self.ny             = obj.read('Grid', 'nr of points y')[0]
        self.nz             = obj.read('Grid', 'nr of points z')[0]
        self.xvector        = obj.read('Grid', 'x-vector')
        self.yvector        = obj.read('Grid', 'y-vector')
        self.zvector        = obj.read('Grid', 'z-vector')

        # Read density
        rho = obj.read('Density', 'KSCED')
        if rho is None: rho = obj.read('Density', 'Active')
        if rho is None: rho = obj.read(type, 'Density')
        if rho is None: rho = obj.read(type, 'Fitdensity')
        self.rho = rho.reshape(self.nz,self.ny,self.nx)
        self.rho = self.rho.transpose()
        rhotot = self.rho.sum()

    @property
    def x(self):
        temp = zeros((self.nx), dtype=float)
        for i in range(self.nx):
            temp[i] = self.start_point[0] + (self.xvector*i)[0]
        return temp

    @property
    def y(self):
        temp = zeros((self.ny), dtype=float)
        for i in range(self.ny):
            temp[i] = self.start_point[1] + (self.yvector*i)[1]
        return temp

    @property
    def z(self):
        temp = zeros((self.nz), dtype=float)
        for i in range(self.nz):
            temp[i] = self.start_point[2] + (self.zvector*i)[2]
        return temp
