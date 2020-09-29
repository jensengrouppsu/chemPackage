from os.path import join as path_join

try:
    from coordinate_functions import calc_bonds, calc_dist
    from coordinate_functions import minmax_pdist, minmax_cdist, translate
    from coordinate_functions import translate, rotate
except ImportError:
    from prep import compile_on_fly
    compile_on_fly(path_join('chem', 'f2py'), 'coordinate_functions')
    from coordinate_functions import calc_bonds, calc_dist
    from coordinate_functions import minmax_pdist, minmax_cdist
    from coordinate_functions import translate, rotate

try:
    from tests import find_equal
except ImportError:
    from prep import compile_on_fly
    compile_on_fly(path_join('chem', 'f2py'), 'tests')
    from tests import find_equal

try:
    from drawingMath import drawing_math
    from drawingMath import drawing_math_vectors
except ImportError:
    from prep import compile_on_fly
    compile_on_fly(path_join('chem', 'f2py'), 'drawingMath')
    from drawingMath import drawing_math
    from drawingMath import drawing_math_vectors
    
__all__ = [
           'calc_bonds',
           'calc_dist',
           'minmax_pdist',
           'minmax_cdist',
           'translate',
           'find_equal',
           'drawing_math'
           'drawing_math_vectors'
          ]
