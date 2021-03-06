from __future__ import print_function, division
import sys, os

'''This is a collection of subroutines that are used at the start of
execution of a program, typically when processing files or command-
line arguements.
'''

def abs_file_path(filename, env=False):
    '''\
    This function takes a *filename* and returns the absolute path.

    The reason this was written is that :py:func:`os.path.abspath`
    can convert a relative path, :py:func:`os.path.expandvars`
    can expand a shell variable, and :py:func:`os.path.expanduser`.
    understands ~, but none of these does all three.  This function
    piggybacks the three of these to guaruntee any path will be
    returned absolutely.

    :argument filename:
        The path of the file that you wish to have the absolute
        path of.
    :type filename: str
    :argument env:
        Replace the base part of the path with $HOME, if that is
        where the path is.  The default is False.
    :type env: bool, optional
    :rtype: str
    '''
    from os.path import abspath, expanduser, expandvars
    absfile = abspath(expandvars(expanduser(filename)))

    # Now replace front with $HOME if requested
    if env:
        if os.environ['HOME'] in absfile:
            i = len(os.environ['HOME']) + 1
            # Assemble absolute path, using the $HOME variable
            absfile = os.path.join('$HOME', absfile[i:])

    return absfile


def locate(pattern, root=None):
    '''\
    Locate a file in the directory tree with a UNIX-style pattern.

    :argument pattern:
        The UNIX-style pattern (uses only ``*``, ? or [ ] syntax) to use
        for searching.
    :type pattern: str
    :keyword root:
        The directory to start searching in.  The default is the
        current directory.
    :type root: str, optional
    :rtype: str
    '''
    import fnmatch
    if root is None: root = os.curdir
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def file_safety_check(filename):
    '''Function to check a file and raise an :py:exc:`IOError`
    if it is not "safe."

    :argument filename:
        The file you wish to check the safety of.
    :type filename: str
    :rtype: None
    :exception:
        :py:exc:`IOError`: Raised when a file is not safe (and abort is
        :py:const:`False`).
    '''
    with open(filename) as f:
        pass
