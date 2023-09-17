from __future__ import print_function, division

def read_file(self):
    '''Reads in the file and store where major sections begin.'''

    # Collect all data into memory for data retention
    with open(self.filename) as fl:
        f = tuple([line.rstrip() for line in fl])

    # For an input file, grab start and end of input block and return
    if self.filetype != 'out':
        return f, { 'INPUT START' : 0, 'INPUT END' : len(f) }

    # Otherwise, read the entire output file

    # Define lines that are accociated with various keys or properties.
    # Since these lines may appear more than once, we define three groups:
    # one where we want the first appearence, one where we want the last
    # appearance, and one where we want every appearence
    # The keys are the lines to find, and the first index of the value is the
    # associated propery, and the second is the number to add the the line
    # number to find the start of where to collect the property.
    #
    # Note that since we are defining the entire line, not just a part,
    # we can search these lines as part of a set which are designed for
    # searching.  If we defined a part of the line, we would have to scan each
    # line, being much slower.
    first = {
            #############
            # Input block
            #############

            # End of the input block
            (' ************* D I S C R E T E   '
              'I N T E R A C T I O N   M O D E L *************'):
                                                             ['INPUT END', -1],
            #############
            # Coordinates
            #############

            # Geometry
            '                                 G E O M E T R Y':
                                                                 ['COORDS', 6],
            ###########
            # Technical
            ###########

            # Start time and date
            '                               Runtime Conditions':
                                                                ['RUNTIME', 3],
            # Timing stats
            '                                Timing Statistics':
                                                                 ['TIMING', 5],
    }
    last = { }
    each = {
            ################################
            # Frequency dependent properties
            ################################

            # The frequency header
            '                    *****  Performing calculation at:  *****':
                                                              ['FREQUENCY', 2],
            # The Polarizabilities
            '                              Polarizability Tensor':
                                                                    ['POL', 0],
            # Induced dipoles and charges
            '         Atom                Dipole                  Charge':
                                                                ['CPIMDIP', 0],
            # Induced dipoles
            '         Atom                Dipole':
                                                                 ['PIMDIP', 0],
            # Efficiencies
            '                         Efficiencies and Cross-Sections':
                                                                    ['EFF', 4],
    }


    #search_lines = set(first.keys()+last.keys()+each.keys())
    search_lines = set(sorted(first)+sorted(last)+sorted(each))
    indices = {}
    for i, line in enumerate(f):
        if line in search_lines:
            # If the line is in the first dict, store then remove from search
            if line in first:
                indices[first[line][0]] = i+first[line][1]
                del first[line]
                search_lines.remove(line)
            # If the line is in a last dict, store and overwrite
            elif line in last:
                indices[last[line][0]] = i+last[line][1]
            # Otherwise, append this value to the index
            else:
                try:
                    indices[each[line][0]].append(i+each[line][1])
                except KeyError:
                    indices[each[line][0]] = []
                    indices[each[line][0]].append(i+each[line][1])

    # Call the input start the beginning of the file
    indices['INPUT START'] = 0
    return f, indices
