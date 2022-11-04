# chemPackage
Python package created to manipulate NWChem and ADF input and output files along with plain XYZ files

installation of the package is done from the master directory through pip by running
    
    pip install .

The chemPackage is dependent on prep and mfunc, to install them run the previous command in both prepContainer and mfuncContainer.


## Dependencies
As of now, the only dependency outside the normal realm of packages like numpy is *natsort*. This package can be found and installed through pip or conda.

## Additional Notes
**This package is still experimental. There are files and routines that still do not have unit tests and have not been verified. If you use any functions that do not have a corresponding unit test, there is a high chance that the results may not be correct. If any issues are found, please open an issue and it will be addressed.**

## example
The main use of the chemPackage all centers around the collect function
    
    In [1]: from chemPackage import collect

    In [2]: f = collect('unitTests/sampleCoords.xyz')

    In [3]: f.coordinates
    Out[3]:
    array([[ 0.  ,  0.  ,  0.  ],
           [ 0.75,  0.75,  0.  ],
           [-0.75,  0.75,  0.  ],
           [ 0.  ,  0.  , -5.  ],
           [ 0.  ,  0.  ,  5.  ],
           [ 0.  , 10.  ,  0.  ],
           [ 0.  ,  0.  ,  7.  ]])
