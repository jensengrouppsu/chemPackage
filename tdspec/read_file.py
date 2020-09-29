from __future__ import print_function, division

def read_file(self):
    '''Reads in the file and store where major sections begin.'''
 
    # Collect all data into memory for data retention
    with open(self.filename) as fl:
        f = tuple([line.rstrip() for line in fl])

    # For an input file, grab start and end of input block and return
    #if self.filetype != 'out':
    #    return f, { 'INPUT START' : 0, 'INPUT END' : len(f) }

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
            ###################
            # Input information
            ###################

            # Start of input information
            ' # TDSPEC OPA Calculation':
                                                              ['OPA START', 0],
            ' # TDSPEC RRS Calculation':
                                                              ['RRS START', 0],
            ' # TDSPEC TPA Calculation':
                                                              ['TPA START', 0],
            ' # TDSPEC RHRS Calculation':
                                                             ['RHRS START', 0],
            ' # TDSPEC DR-SFG Calculation':
                                                           ['DR-SFG START', 0],

            #############
            # End of file
            #############
            ' # TDSPEC OPA Spectrum':
                                                                ['OPA END', 0],
            ' # TDSPEC RRS Spectrum':
                                                                ['RRS END', 0],
            ' # TDSPEC TPA Spectrum':
                                                                ['TPA END', 0],
            ' # TDSPEC RHRS Spectrum':
                                                               ['RHRS END', 0],
            ' # TDSPEC DR-SFG Spectrum':
                                                             ['DR-SFG END', 0],
    }
    last = { }
    each = {
            ########################
            # Calculation Properties
            ########################
            ' Excitation Frequency: ':
                                                      ['EXCITATION ENERGY', 0],
            ##############################################
            # Polarizabilities and normal mode information
            ##############################################
            ' Mode (cm-1)':
                                                            ['NORMAL MODE', 0],
            ' alpha_xx':
                                                         ['POLARIZABILITY', 0],
            ' A_x,xx':
                                                               ['A-TENSOR', 0],
            ' As_x,xx':
                                                              ['As-TENSOR', 0],
            ' G_xx':
                                                               ['G-TENSOR', 0],
            ' Gs_xx':
                                                              ['Gs-TENSOR', 0],
            ' C_xx,xx':
                                                               ['C-TENSOR', 0],
            ' D_xx,x':
                                                               ['D-TENSOR', 0],
            ' Ds_x,xx':
                                                              ['Ds-TENSOR', 0],
            # Zhongwei: Hyperpolarizabilities
            ' beta_xxx':
                                                    ['HYPERPOLARIZABILITY', 0],
    }

    search_lines = set(first.keys()+last.keys()+each.keys())
    indices = {}
    datalines = []
    for i, line in enumerate(f):
        if '#' not in line:
            datalines.append(i)
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
        else:
            # Modes require partial line matches
            for word in each.keys():
                if line[:len(word)] == word:
                    try:
                        indices[each[word][0]].append(i+each[word][1])
                    except KeyError:
                        indices[each[word][0]] = []
                        indices[each[word][0]].append(i+each[word][1])
            

    # Determine where the data starts and ends
    indices['DATA START'] = datalines[0]
    indices['DATA END'] = datalines[-1]

    return f, indices
