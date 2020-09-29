from __future__ import print_function, division
from numpy import array
import os

def collect_input(self, f, indices):
    '''Collect from the input section of the file.'''

    # Keeps track of the index each keyword is found at
    index = {}

    # Create an upper-case and stripped duplicate of the input block.
    # Remove comments (#).
    si, ei = indices['INPUT START'], indices['INPUT END']
    f = list(f)
    search = f[:ei]
    for i, x in enumerate(search[si:], si):
        search[i] = x.upper().lstrip()
        search[i] = search[i].partition('#')[0]
        f[i] = f[i].partition('#')[0]
    f = tuple(f)

    # Let's first look for block-type keywords.  The block-type has the
    # keyword followed by a parameter block terminated with 'END'. 
    # We'll try each known block-type keyword, and if it located find
    # the END keyword and collect.  Otherwise, we'll move to the next
    # keyword.
    en = enumerate
    for keyword in self.blockkeys:
        s = next((i for i, x in en(search[si:], si) if x == keyword), -1)
        if s == -1:
            pass
        else:
            s += 1
            temp = iter(i for i, x in en(search[s:], s) if x == 'END')
            e = next(temp)
#   jbb5516: removing this loop for DIMPAR, does not exist for ADF
  #          if keyword == 'DIMPAR':
  #              try:
  #                  # Old way where we explicitly state the number of elements
  #                  # in the input file.
  #                  for k in range(int(f[s])+1):
  #                      e = next(temp)
  #              except ValueError:
  #                  # New way where we define an element block with the 'element'
  #                  # keyword and implcitly determine how many there are.
  #                  complete = False
  #                  while not complete:
  #                      elements = 0
  #                      ends = 0
  #                      for k in f[s:e]:
  #                          if k.lower().find('element') > -1:
  #                              elements += 1
  #                          if k.lower().find('end') > -1:
  #                              ends += 1
  #                      if elements == ends:
  #                          complete = True
  #                      else:
  #                          e = next(temp)
                    #print(f[e])
                    #print(f[s:e])
            self.key[keyword] = tuple(x for x in f[s:e])
            __determine_subkeys(self, keyword, s, e, search)
            index[keyword] = (s-1, e)
            if keyword == 'DIMPAR': __DIM(self, s, e, search, f)

    # We need to set this because it is possible to do an unrestricted
    # singlet state calculation without specifying mult in the NWChem
    # input file.  This will fix a nasty problem with the collection
    # of excitations...
    self.spin_multiplicity = 1

    # Now, lets locate the line-type keywords. The line-type has the
    # keyword and parameters on the same line.
    for keyword in self.linekeys:
        l = len(keyword)
        ar = [i for i, x in en(search[si:], si) if keyword in x[:l]]
        for ix in ar:
            # Split this line, then add parameters after keyword to key
            line = ' '.join(f[ix].lstrip().split()[1:])
            # Can have multiple instances of the keyword
            try:
                self.key[keyword].append(line)
            except KeyError:
                self.key[keyword] = []
                self.key[keyword].append(line)
            if keyword == 'TITLE':
                self.title = self.key['TITLE'][-1].strip('"')
            elif keyword == 'MULT':
                try:
                    self.spin_multiplicity = int(self.key['MULT'][0])
                except ValueError: #No value given after MULT defaults to 1
                    self.spin_multiplicity = 1
            # Record where found
            index[keyword] = (ix, None)

    # Now, let's find the single-type keywords.  The single-type
    # is just a keyword with no parameters.
    for keyword in self.singlekeys:
        ix = next((i for i, x in en(search[si:], si) if x == keyword), -1)
        if ix == -1:
            pass
        else:
            self.key[keyword] = True
            # Record where found
            index[keyword] = (ix, None)

    # Next, collect the mixed types.  Can be line only or block.
    # Lines are kept in element one, block in element two.
    for keyword in self.mixedkeys:
        ix = next((i for i, x in en(search[si:], si) if keyword in x), -1)
        if ix == -1:
            pass
        else:
            # Make sure that the keyword is the first part of the line
            if search[ix].split()[0] != keyword:
                continue
            self.key[keyword] = ['', tuple()]
            # Add the line part, if any
            self.key[keyword][0]= ' '.join(f[ix].lstrip().split()[1:])
            # Now add the block part.
            s = ix + 1
            # Locate END for this block
            e = next(i for i, x in en(search[s:], s) if x == 'END')
            # Make the parameters a tuple in this key
            self.key[keyword][1] = tuple(x for x in f[s:e])
            # Locate possible subkeys
            __determine_subkeys(self, keyword, s, e, search)
            # GEOMETRY specific handling
            if keyword == 'GEOMETRY':
                from .coordinates import geometry_block
                geometry_block(self, s, e, f)
            # Record where found
            index[keyword] = (ix, e)
            # Make the entry a tuple
            self.key[keyword] = tuple(self.key[keyword])

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
# This seems to not work for nwchem dim, should maybe check later, but 
# for now can simply use external xyz and it works
                self.dim_atoms = array([x.split()[0] for x in f[sx:ex]])
                self.dim_coordinates = array(
                         [x.split()[1:] for x in search[sx:ex]], dtype=float)

def __determine_subkeys(self, k, s, e, search):
    '''Look through the blocks and remember important subkeys.'''

    if k == 'PROPERTY':
        if next((x for x in search[s:e] if 'RESPONSE' in x), False):
            self.subkey.add('RESPONSE')
        if next((x for x in search[s:e] if 'AORESPONSE' in x), False):
            self.subkey.add('AORESPONSE')
        if next((x for x in search[s:e] if 'DAMPING' in x), False):
            self.subkey.add('DAMPING')
    elif k == 'GEOMETRY':
        if next((x for x in search[s:e] if 'SYMMETRY' in x), False):
            self.subkey.add('SYMMETRY')
    if k == 'DFT':
        if next((x for x in search[s:e] if 'ODFT' in x), False):
            self.subkey.add('ODFT')
       #     pass
        if next((x for x in search[s:e] if 'MULT' in x), False):
            # MULT 1 or MULT with no value is still restricted
            try:
                mult = next((x for x in search[s:e] if 'MULT' in x)).split()[1]
                if int(mult) != 1:
                    self.subkey.add('ODFT')
                else:
                    raise IndexError
            except IndexError:
                pass 
    if k == 'TDDFT':
        if next((x for x in search[s:e] if 'CDSPECTRUM' in x), False):
            self.subkey.add('CDSPECTRUM')
        if next((x for x in search[s:e] if 'VELOCITY' in x), False):
            self.subkey.add('VELOCITY')
    if k == 'SCF':
        if next((x for x in search[s:e] if 'UHF' in x), False):
            self.subkey.add('UHF')
        if next((x for x in search[s:e] if 'RHF' in x), False):
            self.subkey.add('RHF')
        if next((x for x in search[s:e] if 'ROHF' in x), False):
            self.subkey.add('ROHF')
    if k == 'TCE':
        if next((x for x in search[s:e] if 'LCCD' in x), False):
            self.subkey.add('LCCD')
        if next((x for x in search[s:e] if 'CCD' in x), False):
            self.subkey.add('CCD')
        if next((x for x in search[s:e] if 'LCCSD' in x), False):
            self.subkey.add('LCCSD')
        if next((x for x in search[s:e] if 'CCSD' in x), False):
            self.subkey.add('CCSD')
        if next((x for x in search[s:e] if 'CCSD_ACT' in x), False):
            self.subkey.add('CCSD_ACT')
        if next((x for x in search[s:e] if 'LR-CCSD' in x), False):
            self.subkey.add('LR-CCSD')
        if next((x for x in search[s:e] if 'CC2' in x), False):
            self.subkey.add('CC2')
        if next((x for x in search[s:e] if 'CCSDT' in x), False):
            self.subkey.add('CCSDT')
        if next((x for x in search[s:e] if 'CCSDTA' in x), False):
            self.subkey.add('CCSDTA')
        if next((x for x in search[s:e] if 'CCSDTQ' in x), False):
            self.subkey.add('CCSDTQ')
        if next((x for x in search[s:e] if 'CCSD(T)' in x), False):
            self.subkey.add('CCSD(T)')
        if next((x for x in search[s:e] if 'CCSD[T]' in x), False):
            self.subkey.add('CCSD[T]')
        if next((x for x in search[s:e] if 'CR-CCSD[T]' in x), False):
            self.subkey.add('CR-CCSD[T]')
        if next((x for x in search[s:e] if 'CR-CCSD(T)' in x), False):
            self.subkey.add('CR-CCSD(T)')
        if next((x for x in search[s:e] if 'CCSD(2)_T' in x), False):
            self.subkey.add('CCSD(2)_T')
        if next((x for x in search[s:e] if 'CCSD(2)_TQ' in x), False):
            self.subkey.add('CCSD(2)_TQ')
        if next((x for x in search[s:e] if 'CCSDT(2)_Q' in x), False):
            self.subkey.add('CCSDT(2)_Q')
        if next((x for x in search[s:e] if 'LR-CCSD(T)' in x), False):
            self.subkey.add('LR-CCSD(T)')
        if next((x for x in search[s:e] if 'LR-CCSD(TQ)-1' in x), False):
            self.subkey.add('LR-CCSD(TQ)-1')
        if next((x for x in search[s:e] if 'CREOMSD(T)' in x), False):
            self.subkey.add('CREOMSD(T)')
        if next((x for x in search[s:e] if 'CREOM(T)AC' in x), False):
            self.subkey.add('CREOM(T)AC')
        if next((x for x in search[s:e] if 'QCISD' in x), False):
            self.subkey.add('QCISD')
        if next((x for x in search[s:e] if 'CISD' in x), False):
            self.subkey.add('CISD')
        if next((x for x in search[s:e] if 'CISDT' in x), False):
            self.subkey.add('CISDT')
        if next((x for x in search[s:e] if 'CISDTQ' in x), False):
            self.subkey.add('CISDTQ')
        if next((x for x in search[s:e] if 'MBPT2' in x), False):
            self.subkey.add('MP2')
        if next((x for x in search[s:e] if 'MP2' in x), False):
            self.subkey.add('MP2')
        if next((x for x in search[s:e] if 'MBPT3' in x), False):
            self.subkey.add('MP3')
        if next((x for x in search[s:e] if 'MP3' in x), False):
            self.subkey.add('MP3')
        if next((x for x in search[s:e] if 'MBPT4' in x), False):
            self.subkey.add('MP4')
        if next((x for x in search[s:e] if 'MP4' in x), False):
            self.subkey.add('MP4')
