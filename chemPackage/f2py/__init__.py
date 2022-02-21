from os.path import join as path_join

try:
    from .coordinate_functions import calc_bonds, calc_dist
    from .coordinate_functions import minmax_pdist, minmax_cdist, translate
    from .coordinate_functions import translate, rotate
except ImportError:
    from prep import compile_on_fly
    compile_on_fly(path_join('chemPackage', 'f2py'), 'coordinate_functions')
    from .coordinate_functions import calc_bonds, calc_dist
    from .coordinate_functions import minmax_pdist, minmax_cdist
    from .coordinate_functions import translate, rotate

try:
    from .tests import find_equal
except ImportError:
    from prep import compile_on_fly
    compile_on_fly(path_join('chemPackage', 'f2py'), 'tests')
    from .tests import find_equal

try:
    from .drawingMath import drawing_math
    from .drawingMath import drawing_math_ret
    from .drawingMath import drawing_math_ret_magnetic
    from .drawingMath import drawing_math_vectors
    from .drawingMath import drawing_math_vectors_ret
    from .drawingMath import drawing_math_vectors_ret_magnetic
except ImportError:
    from prep import compile_on_fly
    compile_on_fly(path_join('chemPackage', 'f2py'), 'drawingMath')
    from .drawingMath import drawing_math
    from .drawingMath import drawing_math_ret
    from .drawingMath import drawing_math_ret_magnetic
    from .drawingMath import drawing_math_vectors
    from .drawingMath import drawing_math_vectors_ret
    from .drawingMath import drawing_math_vectors_ret_magnetic
   
__all__ = [
           'calc_bonds',
           'calc_dist',
           'minmax_pdist',
           'minmax_cdist',
           'translate',
           'find_equal',
           'drawing_math',
           'drawing_math_ret',
           'drawing_math_ret_magnetic',
           'drawing_math_vectors',
           'drawing_math_vectors_ret'
           'drawing_math_vectors_ret_magnetic'
          ]
