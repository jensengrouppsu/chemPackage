from __future__ import print_function, division
from datetime import datetime, timedelta
from os.path import basename
from numpy import where

def collect_technical(self, f, indices):
    '''Collect termination status, timing data, and computation location.'''

    ##############################
    # DETERMINE TERMINATION STATUS
    ##############################

    # The dermination status appears 5 lines before 'Timing Statistics'.
    self.termination = 'Unknown Calculation Failure'
    if 'ERROR' in indices:
        self.termination = 'ERROR: ' + f[indices['ERROR']].strip()
    elif 'TIMING' in indices:
        self.termination = f[indices['TIMING']-5].strip()

    #####################
    # COLLECT TIMING DATA
    #####################

    # Collect the date and time the calculation started at.
    # Use only first instance this is found
    try:
        tl = next((x for x in f[indices['INPUT END']:] if 'RunTime:' in x), None)
    except KeyError:
        tl = None
    if tl is not None:
        # Grab date and time from MonDD-YYYY HH:MM:SS format
        # and place into datetime object.
        if tl.split()[1] == 'development':
            tp = tl.split()[4] + ' ' + tl.split()[5]
        elif tl.split()[1] == 'RunTime:':
            tp = tl.split()[2] + ' ' + tl.split()[3]
        else:
            tp = tl.split()[3] + ' ' + tl.split()[4]
        self.start = datetime.strptime(tp, '%b%d-%Y %H:%M:%S')
        
    #########################################
    # DETERMINE WHERE CALCULATION RAN AND HOW
    #########################################

    # Parallel info tells the hostname of the computer and # procs.
    if 'OLD PARALLELIZATION' in indices:
        # The last number three above header is the # procs
        self.nprocs = int(f[indices['OLD PARALLELIZATION']-3].split()[-1])
        # Last word two above header is the hostname
        self.host = f[indices['OLD PARALLELIZATION']-2].split()[-1]

    elif 'NEW PARALLELIZATION' in indices:
        # Look only at last entry for id.  This will give one less
        # than the number of nodes as well.  This is the line before
        # The terminating '=========='
        s = indices['NEW PARALLELIZATION']
        ix = next(i for i, x in enumerate(f[s:], s) if '=======' in x)
        tp = f[ix-1]
        # Second word is host
        self.host = tp.split()[1]
        # First number is rank of this process
        self.nprocs = int(tp.split()[0]) + 1

    # Now make the computer names uniform using predefined names
    if self.host is not None:
        for name in self._computer_names:
            if name.lower() in self.host: self.host = name
