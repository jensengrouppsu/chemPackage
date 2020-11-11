#! /usr/bin/env python

from __future__ import print_function

def open_file(prefix):
    '''Opens a VTR file for writting and returns the related file object'''

    filename = prefix+'.vtr'
    fd = open(filename, 'w')
    print ('<VTKFile type="RectilinearGrid" version="0.1" format="ascii">',
           file=fd)
    return fd

def write_mesh_3D(fd, x, y, z):
    '''Writes a 3D mesh to file object fd'''

    print (' <RectilinearGrid WholeExtent="1 {0} 1 {1} 1 {2}">'.format(
           len(x), len(y), len(z)), file=fd)
    print ('  <Piece Extent="1 {0} 1 {1} 1 {2}">'.format(
           len(x), len(y), len(z)), file=fd)
    print ('   <Coordinates>', file=fd)
    print ('    <DataArray type="Float32" Name="X_COORDINATES" NumberOfComponents="1" format="ascii">',
           file=fd)
    for i in range(len(x)):
        if (i+1 == len(x)) or ((i+1)%5 == 0):
            end = '\n'
        else:
            end = '  '
        print ('{0:13.5f}'.format(x[i]), end=end, file=fd)
    print ('    </DataArray>', file=fd)
    print ('    <DataArray type="Float32" Name="Y_COORDINATES" NumberOfComponents="1" format="ascii">',
           file=fd)
    for i in range(len(y)):
        if (i+1 == len(y)) or ((i+1)%5 == 0):
            end = '\n'
        else:
            end = '  '
        print ('{0:13.5f}'.format(y[i]), end=end, file=fd)
    print ('    </DataArray>', file=fd)
    print ('    <DataArray type="Float32" Name="Z_COORDINATES" NumberOfComponents="1" format="ascii">',
           file=fd)
    for i in range(len(z)):
        if (i+1 == len(z)) or ((i+1)%5 == 0):
            end = '\n'
        else:
            end = '  '
        print ('{0:13.5f}'.format(z[i]), end=end, file=fd)
    print ('    </DataArray>', file=fd)
    print ('   </Coordinates>', file=fd)
    print ('   <PointData>', file=fd)

def write_scalar_3D(fd, name, field):
    '''Writes a scalar onto a 3D grid'''

    print ('         <DataArray type="Float32" Name="'
          +name+'" NumberOfComponents="1" format="ascii">', file=fd)
    n = 0
    nmax = len(field[0]) * len(field[1]) * len(field[0][0])
    for i in range(len(field[0][0])):
        for j in range(len(field[0])):
            for k in range(len(field)):
                n += 1
                if (n == nmax) or (n%5 == 0):
                    end = '\n'
                else:
                    end = '  '
                print ('{0:13.5E}'.format(field[k][j][i]), end=end, file=fd)
    print ('         </DataArray>', file=fd)

def close_file(fd):
    '''Closes the vtr file object'''

    print ('     </PointData>', file=fd)
    print ('   </Piece>', file=fd)
    print ('  </RectilinearGrid>', file=fd)
    print ('</VTKFile>', file=fd)
    fd.close()

def make_3D_vtr_file(prefix, x, y, z, field):
    '''Performs necessary routines to create a 3D vtr file of a scalar'''

    fd = open_file(prefix)
    write_mesh_3D(fd, x, y, z)
    write_scalar_3D(fd, prefix, field)
    close_file(fd)
