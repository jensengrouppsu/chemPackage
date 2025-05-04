import numpy as np

# Modified from adf/input_block.py  --Gaohe
def collect_input(self, f, indices):
    '''Collect from the input section of the file.'''

    # Keeps track of the index each keyword is found at
    index = {}

    # Create an upper-case and stripped duplicate of the input block.
    # Remove comments (:: or !).
    si, ei = indices['INPUT START'], indices['INPUT END']
    f = list(f)
    search = f[:ei]
    for i, x in enumerate(search[si:], si):
        search[i] = x.upper().lstrip()
        search[i] = search[i].partition('::')[0]
        search[i] = search[i].partition('!')[0]
        f[i] = f[i].partition('::')[0]
        f[i] = f[i].partition('!')[0]
    f = tuple(f)

    # Let's first look for block-type keywords.  The block-type has the
    # keyword followed by a parameter block terminated with 'END'. 
    # We'll try each known block-type keyword, and if it located find
    # the END keyword and collect.  Otherwise, we'll move to the next
    # keyword.
    for keyword in self.blockkeys:
        s = next((i for i, x in enumerate(search[si:], si) if x == keyword), -1)
        if s == -1:
            pass
        else:
            s += 1
            # Locate END for this block
            e = next(i for i, x in enumerate(search[s:], s) if x == 'END')
            # Make the parameters a tuple in this key
            self.key[keyword] = tuple(x for x in f[s:e])
            # Locate possible subkeys
            __determine_subkeys(self, keyword, s, e, search)
            # Record where found
            index[keyword] = (s-1, e)
            # DIMQMPAR specific handling
            if keyword == 'DIMPAR': __DIM(self, s, e, search, f)

    # Now, lets locate the line-type keywords. The line-type has the
    # keyword and parameters on the same line.
    for keyword in self.linekeys:
        l = len(keyword)
        ar = [i for i, x in enumerate(search[si:], si) if keyword == x[:l]]
        for ix in ar:
            # Split this line, then add parameters after keyword to key
            line = ' '.join(f[ix].lstrip().split()[1:])
            # Multiple instances possible
            try:
                self.key[keyword].append(line)
            except KeyError:
                self.key[keyword] = []
                self.key[keyword].append(line)
            if keyword == 'TITLE': self.title = self.key['TITLE']
            # Record where found
            index[keyword] = (ix, None)

    # Now, let's find the single-type keywords.  The single-type
    # is just a keyword with no parameters.
    for keyword in self.singlekeys:
        ix = next((i for i, x in enumerate(search[si:], si) if x == keyword), -1)
        if ix == -1:
            pass
        else:
            self.key[keyword] = True
            # Record where found
            index[keyword] = (ix, None)

    # Next, collect the mixed types.  Can be line only or block.
    # Lines are kept in element one, block in element two.
    # FIXME: Some of these keys have been moved to engine settings, commented out for now.
    # for keyword in self.mixedkeys:
    #     ix = next((i for i, x in enumerate(search[si:], si) if keyword in x), -1)
    #     if ix == -1:
    #         pass
    #     else:
    #         # Skip EFIELD if this word is in the DIMQM block
    #         if keyword == 'EFIELD' and 'DIMQM' in self.key:
    #             if index['DIMQM'][0] < ix < index['DIMQM'][1]: continue
    #         # Make a two element list
    #         self.key[keyword] = ['', tuple()]
    #         # Add the line part, if any
    #         self.key[keyword][0] = ' '.join(f[ix].lstrip().split()[1:])
    #         # Now add the block part if it exists.  Not always true
    #         # For EFIELD or GEOMETRY
    #         bl = '&' in self.key[keyword][0] or not self.key[keyword][0]
    #         if ((keyword in ('GEOMETRY', 'EFIELD') and bl) or
    #             (keyword in ('FRAGMENTS', 'ATOMS'))):
    #             s = ix + 1
    #             # Locate END for this block
    #             e = next(i for i, x in enumerate(search[s:], s) if x == 'END')
    #             # Make the parameters a tuple in this key
    #             self.key[keyword][1] = tuple(x for x in f[s:e])
    #             # Locate possible subkeys
    #             __determine_subkeys(self, keyword, s, e, search)
    #             # ATOMS specific handling
    #             if keyword == 'ATOMS':
    #                 from .coordinates import atom_block
    #                 atom_block(self, s, e, search, f, indices)
    #             # Record where found
    #             index[keyword] = (ix, e)
    #         else:
    #             # Record where found
    #             index[keyword] = (ix, None)
    #         # Make the entry a tuple
    #         self.key[keyword] = tuple(self.key[keyword])

    # Last, keep a record of the order of the inputs
    self._input_order = tuple(sorted(index, key=index.get))


def __DIM(self, s, e, search, f):
    '''Support for DIM-specific things.'''
    # Import DIM coordinates from an external xyz file or block.
    # Only do if not an output file, since we will collect later otherwise.
    if self.filetype != 'out':
        # First see if the XYZ block is embedded in the file
        sx = next((i for i, x in enumerate(search[s:], s) if x == 'XYZ'), -1)
        if sx == -1:
            # If it is not embedded in the file, see if it is in another file
            try:
                xyzfile = next(x.lstrip() for x in f[s:e] if '.xyz' in x)
            except StopIteration:
                # Last, check if it is in a GROUP block
                pass
            # Read coordinates from another file
            else:
                from .. import collect
                try:
                    xyz = collect(xyzfile)
                except IOError:
                    pass
                else:
                    self.ndim = xyz.natoms
                    self.dim_atoms = xyz.atoms
                    self.dim_coordinates = xyz.coordinates
        # Read coordinates from an XYZ block
        else:
            ex = next(i for i, x in enumerate(search[sx:],sx) if x == 'SUBEND')
            if '.xyz' in f[sx+1]:
                xyzfile = f[sx+1].lstrip()
                from .. import collect
                try:
                    xyz = collect(xyzfile)
                except IOError:
                    pass
                else:
                    self.ndim = xyz.natoms
                    self.dim_atoms = xyz.atoms
                    self.dim_coordinates = xyz.coordinates
            else:
                sx = sx + 2
                self.dim_atoms = np.array([x.split()[0] for x in f[sx:ex]])
                self.dim_coordinates = np.array(
                         [x.split()[1:] for x in search[sx:ex]], dtype=float)

def __determine_subkeys(self, k, s, e, search):
    '''Look through the blocks and remember important subkeys.'''

    if k == 'GEOMETRY':
        if self.key[k][0] and 'FREQUENCIES' in self.key[k][0]:
            self.subkey.add('FREQUENCIES')
        elif next((x for x in search[s:e] if 'FREQUENCIES' in x), False):
            self.subkey.add('FREQUENCIES')
    # elif k == 'ANALYTICALFREQ': self.subkey.add('FREQUENCIES')
    elif k == 'AORESPONSE':
        if next((x for x in search[s:e] if 'LIFETIME' in x), False):
            self.subkey.add('LIFETIME')
        if next((x for x in search[s:e] if 'OPTICALROTATION' in x), False):
            self.subkey.add('OPTICAL ROTATION')
        if next((x for x in search[s:e] if 'RAMAN' in x), False):
            self.subkey.add('RAMAN')
        if next((x for x in search[s:e] if 'VROA' in x), False):
            self.subkey.add('VROA')
        if next((x for x in search[s:e] if 'QUADRUPOLE' in x), False):
            self.subkey.add('A-TENSOR')
        if next((x for x in search[s:e] if 'TWONPLUSONE' in x or 'BETA' in x), False):
            self.subkey.add('TWONPLUSONE' or 'BETA')
        #if next((x for x in search[s:e] if 'BETA'), False):
        #    self.subkey.add('BETA')
        if next((x for x in search[s:e] if 'GAMMA' in x), False):
            self.subkey.add('GAMMA')
        # TODO: I didn't find any places that use this subkey, and it looks wrong
        # if next((x for x in search[s:e] if 'FREQRANGE' in x), False):
        #     self.subkey.add('FREQRANGE')
    elif k == 'RESPONSE':
        if next((x for x in search[s:e] if 'RAMAN' in x), False):
            self.subkey.add('RAMAN')
        if next((x for x in search[s:e] if 'HYPERPOL' in x), False):
            self.subkey.add('HYPERPOL')
    elif k == 'EXCITATION':
        if next((x for x in search[s:e] if 'EXACT' in x), False):
            self.subkey.add('EXACT')
        if next((x for x in search[s:e] if 'TD-DFTB' in x), False):
            self.subkey.add('TD-DFT+TB')
    elif k == 'DIMQM':
        if next((x for x in search[s:e] if 'CPIM' in x), False):
            self.subkey.add('CPIM')
    elif k =='FDE':
        if next((x for x in search[s:e] if 'FREEZEANDTHAWCYCLES' in x), False):
            self.subkey.add('FDE-RELAX')
        if next((x for x in search[s:e] if 'RELAXCYCLES' in x), False):
            self.subkey.add('FDE-RELAX')
