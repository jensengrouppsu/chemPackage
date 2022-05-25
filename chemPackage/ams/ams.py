from ..chemdata import ChemData
from ..errorclass import CollectionError
from numpy import array, zeros, transpose


# Create the AMS class to grab and hold all the information in the file
class AMS(ChemData):
    '''Read and contain data from AMS files
    
    This class is used to read all information from either an AMS input or
    output file.  It stores this data in formats suitable for processing.  

    The only required argument for collection is the name of the file.
 
    The data available is listed below.  If the data did not exist in the
    AMS file, then the variable will default to None (or an empty set for
    subkey and calctype, or empty dict for key).
    
    '''
    # Instantiate the AMS class
    def __init__(self, name):
        '''Set up the AMS class.  The only argument is the filename.'''
        import os

        ChemData.__init__(self)
        self.program = 'AMS'
        self.project='all'
    
        # Find extention
        ftype = os.path.splitext(name)[1]
        if ftype not in ('.out', '.run', '.inp'):
            raise ValueError (ftype+' not a recognized extention')
        self.filetype = ftype[1:]
        self.filename = name

        # These keys may or may not be good
        # Create tuples of general keywords
        self.mixedkeys = ('ATOMS', 'EFIELD', 'FRAGMENTS', 'GEOMETRY',
                          'INTEGRATION')
        self.blockkeys = ('ANALYTICALFREQ', 'AORESPONSE', 'BASIS',
                          'CONSTRAINTS', 'DIMQM', 'DIMPAR', 'EXCITATION',
                          'EXCITEDGO', 'EXTERNALS', 'FDE', 'GEOVAR',
                          'GUIBONDS', 'RESPONSE', 'REMOVEFRAGORBITALS', 'SCF',
                          'SOLVATION', 'UNITS', 'VIBRON', 'XC', 'SUBEXCI', 'SUBRESPONSE')
        self.linekeys = ('A1FIT', 'BONDORDER', 'CHARGE', 'CREATE',
                         'NOPRINT', 'PRINT', 'RELATIVISTIC',
                         'SAVE', 'SCANFREQ', 'SYMMETRY', 'THERMO', 'TITLE',
                         'DEPENDENCY')
        self.singlekeys = ('ALLPOINTS', 'BADER', 'FORCEALDA', 'NEWDIIS',
                           'UNRESTRICTED', 'DIFFUSE', 'EXACTDENSITY',
                           'IGNOREOVERLAP', 'AOMAT2FILE', 'STOFIT', 'TOTALENERGY',
                           'AOMAT2FILE', 'ALLOWPARTIALSUPERFRAGS')

    def _collect(self, abort=False):
        '''Collect the AMS data from file.  'abort' will cause collection
        to stop on an error, instead of ignoring the error. '''

        # Save the abort status
        self._abort = abort

        # Conventions used in _collect and helper functions:
        # fi = input block
        # fo = output information block
        # s = start of current collection block
        # e = end of current collection block
        # tp = general temporary storage
        # tl = temporary line storage
        # ln = a split line
        # ix = temporary index array
        # ar = temporary array
        # bl = a temporary boolean
        # sl = a line to search for
        # fn = a temporary function
        # m = the results of a regex search
        if self.project=='all':
            # Read in file
            from .read_file import read_file
            f, indices = read_file(self)
    

            if self.filetype != 'out': return
            
            if 'INITIAL GEOMETRY' in indices:
                from .coordinates import collect_geometry
                collect_geometry(self, f, indices)
    
            if 'OPTIMIZED GEOMETRY' in indices:
                from .coordinates import collect_optimized_geometry
                collect_optimized_geometry(self, f, indices)

            # Collect frequencies and IR intensities
            if "NORMAL MODES" in indices:
                from .vibrations import collect_frequencies
                collect_frequencies(self, f, indices)
                self.calctype.add('FREQUENCIES')

            if "MBH" in indices:
                from .vibrations import collect_frequencies_mbh
                collect_frequencies_mbh(self, f, indices)
                self.calctype.add('FREQUENCIES')

            # Collect polarizability
            if 'AORESPONSE' in indices:
                from .polarizability import collect_polarizability
                collect_polarizability(self, f, indices)

            # Collect G-Tensor
            if 'OPTICAL ROTATION' in indices:
                from .polarizability import collect_opticalrotation
                collect_opticalrotation(self, f, indices)
            
            # Collect A-Tensor
            if 'A-TENSOR' in indices:
                from .polarizability import collect_atensor
                collect_atensor(self, f, indices)

            # Collect C-Tensor 
            if 'C-TENSOR' in indices:
                from .polarizability import collect_ctensor
                collect_ctensor(self, f, indices)

            # Collecct magnetizability
            if 'MAGNETIZABILITY' in indices:
                from .polarizability import collect_magnetizability
                collect_magnetizability(self, f, indices)

            # collect linear response
            if "LINEAR RESPONSE" in indices:
                from .polarizability import collect_linearresponse
                collect_linearresponse(self, f, indices)


