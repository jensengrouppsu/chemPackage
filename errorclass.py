'''A collection of error classes.'''
from .errorclassinfo import __doc__

class CollectionError(Exception):
    '''Error class for file collections.

    :param msg: the message to give to the user.

    Printing off this class as a string will display a pre-made
    error message:

    .. code-block:: python
    
        try:
            filedata = collector(filename)
        except CollectionError as c:
            sys.exit(str(c))

    This keeps implementation simple.  
    
    '''
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg

class DrawError(Exception):
    '''Error class for drawing.'''
    def __init__(self,msg):
        self.msg = msg
    def __str__(self):
        return self.msg

class ChemDataError(Exception):
    '''Error class for ChemData.'''
    def __init__(self, msg):
        self.msg = msg
    def __str__(self):
        return self.msg
