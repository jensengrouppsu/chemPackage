from __future__ import print_function, division
from .errorclass import ChemDataError
from sys import stdout

class InputFiles(object):
    'Extends ChemData class with methods to create and compare input files.'

    def copy_template(self, template, file=None, a1=1, a2=None, charge=None,
                      basis=None):
        '''Prints a new input file based on a template.  Writes to file.'''
        from . import collect
        from os.path import splitext

        # Collect the template just to see the program type
        t = collect(template)

        # Read template as an array
        with open(template) as fl:
            tmplt = [x.rstrip() for x in fl]
        cap = [x.upper() for x in tmplt]

        # Look for the block of coordinates
        if t.program == 'ADF':
            try:
                s = next(i for i, x in enumerate(cap) if 'ATOMS' in x) + 1
            except StopIteration:
                raise ChemDataError ('No ATOMS block in template file.')
            try:
                e = next(i for i, x in enumerate(cap[s:], s) if 'END' in x)
            except StopIteration:
                raise ChemDataError ('No END for ATOMS in template file.')
        elif t.program == 'AMS':
            try:
                s = next(i for i, x in enumerate(cap) if 'ATOMS' in x) + 1
            except StopIteration:
                raise ChemDataError ('No ATOMS block in template file.')
            try:
                e = next(i for i, x in enumerate(cap[s:], s) if 'END' in x)
            except StopIteration:
                raise ChemDataError ('No END for ATOMS in template file.')
        elif t.program == 'NWChem':
            try:
                s = next(i for i, x in enumerate(cap) if 'GEOMETRY' in x) + 1
            except StopIteration:
                raise ChemDataError ('No GEOMETRY block in template file.')
            try:
                e = next(i for i, x in enumerate(cap[s:], s) if 'END' in x)
            except StopIteration:
                raise ChemDataError ('No END for GEOMETRY in template file.')
            ix = next((i for i,x in enumerate(cap[s:e], s) if 'SYMMETRY' in x),
                      None)
            try:
                sym = tmplt[ix]
            except TypeError:
                sym = None
        elif t.program == 'Dalton':
            fullbasis = False
            # The top is either defined by the BASIS or ATOMBASIS keyword.
            try:
                s = next(i for i, x in enumerate(cap) if x == 'BASIS') + 5
                fullbasis = True
            except StopIteration:
                s = next(i for i, x in enumerate(cap) if x == 'ATOMBASIS') + 4
            except StopIteration:
                raise ChemDataError ('No BASIS or ATOMBASIS in template file.')
            try:
                e = next(i for i, x in enumerate(cap[s:], s) if 'DALTON INPUT' in x) - 2
            except StopIteration:
                raise ChemDataError ('No end for geometry block in template file.')

        print('s/e',s,e)

        # Make title format string here, in case it is needed later.
        if t.program == 'ADF':
            tstr = 'TITLE {0}'
        if t.program == 'AMS':
            tstr = ''
        elif t.program == 'NWChem':
            tstr = 'title "{0}"'
        elif t.program == 'Dalton':
            tstr = '{0}'

        # Look for title line for ADF and NWChem
        try:
            ititle = next(i for i, x in enumerate(cap) if 'TITLE' in x)
        except StopIteration:
            ititle = None
            # Reset the start index so that we don't put title inside coords
            s = s - 1
            coordhead = tmplt[s]
        # Replace old title with new one, if it exists
        else:
            # If there is no title to replace, leave blank
            if self.title is None:
                tmplt[ititle] = tstr.format('')
            else:
                tmplt[ititle] = tstr.format(self.title)

        # For Dalton, the second/third line of the input is always the title.
        # self is the coordinate file, t is the template.
        if t.program == 'Dalton':
            if fullbasis:
                ititle = 2
            else:
                ititle = 1
            if self.title is None:
                tmplt[ititle] = tstr.format('')
            else:
                tmplt[ititle] = tstr.format(self.title)

        # If file is None, use standard out
        if file is None:
            out = stdout
            closebool = False
        # Otherwise, try to open the file
        else:
            try:
                out = open(file, 'w')
            # If it fails, then it was already an open file
            except TypeError:
                out = file
                closebool = False
            # If it succeeds, then remember that we must close the file
            else:
                closebool = True

        # Print out the template
        # Everything before the coordinates
        for x in tmplt[0:s]:
            print(x, file=out)

        # If self has a title, but the template did not, add here
        if t.program != 'AMS':
            if ititle is None:
                print(tstr.format(self.title), end='\n\n', file=out)
                print(coordhead, file=out)
        else:
            if ititle is None:
                print(coordhead, file=out)
            

        # The coordinates
        if t.program == 'ADF':
            self.printCoords(mode='num', file=out, a1=a1, a2=a2)
        if t.program == 'AMS':
            self.printCoords(mode=None, file=out, a1=a1, a2=a2)
        elif t.program == 'NWChem':
            self.printCoords(file=out, a1=a1, a2=a2)
            if sym is not None: print(sym, file=out)
        elif t.program == 'Dalton':
            # Make the coordinate block heading
            latom = str(len(str(self.nelements)))
            lcharge = str(len(charge))
            fmt = '{0}={1:<'+latom+'} {2}={3:<'+lcharge+'} {4} {5}'
            print(fmt.format('Atomtypes', self.nelements, 'Charge',
                             charge, 'Angstrom', 'Nosymmetry'), end='\n',
                             file=out)
            self.printCoords(dalton=True, file=out, a1=a1, a2=a2, abasis=basis)

        # Everything after the coordinates
        for x in tmplt[e:]:
            print(x, file=out)

        # Close the file
        if closebool: out.close()


    def inputs_match(self, other, ignore=set()):
        '''See if two input files match, possibly ignoring keys in ignore.

        Ignore is converted to a set internally, so it need not be given
        as one.
        '''
        from difflib import SequenceMatcher
        from pprint import pformat
        from . import collect
        from .chemdata import ChemData
        from .errorclass import CollectionError

        # If other is not a ChemData object, assume it is a file to collect
        if not isinstance(other, ChemData):
            try:
                other = collect(other, raise_err=True)
            # If it doesn't exist or can't be collected, it doesn't match
            except (IOError, CollectionError):
                print("B DOESN'T EXIST")
                return False

        # Populate keys with all key types.
        keys = self.mixedkeys+self.blockkeys+self.linekeys+self.singlekeys

        # Convert ignore and keys into a set.  Remove keys in the ignore set
        keys = set(keys)
        ignore = set(ignore)
        keys.difference_update(ignore)

        # Now make sure that each key is the same.  Uses pformat to change
        # whatever datatype into a string, and SequenctMatcher to compare
        # those strings.
        for key in keys:
            try:
                a = pformat(self.key[key])
            except KeyError:
                a_has_key = False
            else:
                a_has_key = True
            try:
                b = pformat(other.key[key])
            # If one has the key and the other doesn't they don't match
            except KeyError:
                if a_has_key:
                    print("B MISSING KEY", key)
                    return False
                else:
                    continue
            else:
                if not a_has_key: 
                    print("A MISSING KEY", key)
                    return False
            comp = SequenceMatcher(None, a, b)
            # ratio() method returns 1.0 if matches are exact
            if comp.ratio() != 1.0: return False

        # Otherwise, they're the same
        return True
