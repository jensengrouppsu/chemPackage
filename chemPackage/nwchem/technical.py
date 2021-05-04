from __future__ import print_function, division
from datetime import datetime, timedelta
from os.path import basename
from numpy import where

def collect_technical(self, f, indices):
    '''Collect termination status, timing data, and computation location.'''

    ##############################
    # DETERMINE TERMINATION STATUS
    ##############################

    # The termination status appears 5 lines before 'Timing Statistics'.
    err = 'Unknown Calculation Failure'
    if 'TIMING' in indices:
        err = 'NORMAL TERMINATION'
    else:
        # To be implemented sometime in the future
        pass
    self.termination = err

    #######################################################
    # DETERMINE WHERE CALCULATION RAN, NPROC AND START TIME
    #######################################################

    # Find the job info block
    if 'JOB INFO' in indices:
        s = indices['JOB INFO']
        e = indices['MEM INFO']

        # Determine hostname
        self.host = next(x.split()[-1] for x in f[s:e] if 'hostname' in x)
        # Make the computer name uniform using predefined names
        for name in self._computer_names:
            if name.lower() in self.host: self.host = name

        # Determine the number of processors
        self.nprocs = next(i.split()[-1] for i in f[s:e] if 'nproc' in i)
        self.nprocs = int(self.nprocs)

        # Determine the starting time
        self.start = next(x.split()[2:] for x in f[s:e] if 'date' in x)
        self.start = ' '.join(self.start)
        self.start = datetime.strptime(self.start, '%a %b %d %H:%M:%S %Y')
