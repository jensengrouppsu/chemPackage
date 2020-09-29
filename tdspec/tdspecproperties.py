from __future__ import print_function, division
from numpy import array, where, zeros, append, argsort

def collect_tdspec(self, f, indices):
    '''Collect all TDSPEC properties.'''

    #####################
    # Collect TDSPEC data
    #####################
    s = indices['DATA START']
    e = indices['DATA END']

    # Determine what type of calculation was done
    absorption = ('OPA', 'TPA',)
    vibrational = ('RRS', 'RHRS', 'DR-SFG',)
    
    calcabs = False
    calcvib = False

    for ctype in self.calctype:
        if ctype in absorption:
            calcabs = True
        elif ctype in vibrational:
            calcvib = True

    # Absorption type processes
    if calcabs:
        self.tdspec_wavelength = array([x.split()[0] for x in f[s:e]],
                                        dtype=float)
        if 'OPA' in self.calctype:
            self.tdspec_opa = array([x.split()[1] for x in f[s:e]],
                                    dtype=float)
        elif 'TPA' in self.calctype:
            self.tdspec_tpa = array([x.split()[1] for x in f[s:e]],
                                    dtype=float)
    # Vibrational spectroscopies
    elif calcvib:
        self.tdspec_wavenumber = array([x.split()[0] for x in f[s:e]],
                                        dtype=float)
        if 'RRS' in self.calctype:
            self.tdspec_rrs = array([x.split()[1] for x in f[s:e]],
                                    dtype=float)
        elif 'RHRS' in self.calctype:
            self.tdspec_rhrs = array([x.split()[1] for x in f[s:e]],
                                     dtype=float)
        elif 'DR-SFG' in self.calctype:
            self.tdspec_drsfg = array([x.split()[1] for x in f[s:e]],
                                      dtype=float)

def collect_e_frequency(self, f, indices):
    '''Collects the Excitation frequency.'''

    from chem.constants import WAVENUM2HART as W2H
    ar = indices['EXCITATION ENERGY']
    for ix in ar:
        temp = array([W2H(float(f[ix].split()[2]))], dtype=float)
    self.e_frequencies = zeros((len(self.v_frequencies)), dtype=float) + temp
    # ---------------------------------------------------------------------------
    # Zhongwei: two frequencies for HyperRaman
    self.b_e_frequencies = zeros((len(self.v_frequencies)), dtype=float) + temp/2
    # ---------------------------------------------------------------------------

#def collect_timing(self, f, indices):
#    '''Collect the timing info.'''
#    from datetime import datetime, timedelta
#    from mfunc import intround
#
#    # Collect the starting time
#    if 'RUNTIME' in indices:
#        ix = indices['RUNTIME']
#        tp = ' '.join(f[ix].strip().split()[3:]) # The date
#        tp += ' ' + f[ix+1].strip().split()[3]     # Add time to date
#        self.start = datetime.strptime(tp, '%A, %b %d, %Y %X')
#    else:
#        self._raise_or_pass('Could not find DIM runtime conditions')
#
#    # Collect the timings
#    if 'TIMING' in indices:
#        ix = indices['TIMING']
#        self.real_time = timedelta(seconds=intround(f[ix].split()[0]))
#        self.cpu_time  = timedelta(seconds=intround(f[ix].split()[1]))
#        self.routine_times = {}
#        for line in f[ix+1:]:
#            if not line.strip(): continue
#            ln = line.split()
#            # If the routine is multiple words, join the words
#            tp = ' '.join(ln[3:])
#            # 0: Real time, 1: CPU time, 2: # Calls
#            if tp not in self.routine_times:
#                self.routine_times[tp] = (timedelta(seconds=intround(ln[0])),
#                                          timedelta(seconds=intround(ln[1])),
#                                          int(ln[2]))
#            else:
#                # Add the new time to the old time
#                self.routine_times[tp] = list(self.routine_times[tp])
#                self.routine_times[tp][0] += timedelta(seconds=intround(ln[0]))
#                self.routine_times[tp][1] += timedelta(seconds=intround(ln[1]))
#                self.routine_times[tp][2] += int(ln[2])
#                self.routine_times[tp] = tuple(self.routine_times[tp])
#
#def collect_technical(self, f, indices):
#    '''Collect technical info, such as where the job was run and termination'''
#
#    # Look above the timing for an error message
#    if 'TIMING' in indices:
#        ix = indices['TIMING']
#        try:
#            # Look up to 50 lines before the timing for an error
#            ix = next(i for i in xrange(ix, ix-50, -1) if 'ERROR' in f[i])
#        except StopIteration:
#            # If no error was found then we terminated correctly
#            self.termination = 'NORMAL TERMINATION'
#        else:
#            # An error was found, save it
#            self.termination = f[ix].strip().replace('ERROR: ', '')
#        
#    # Find the number of processor
#    if 'RUNTIME' in indices:
#        ix = indices['RUNTIME']
#        
#        # Get the host name
#        self.host = f[ix+2].split()[3]
#        # Get the number of processors
#        self.nprocs = int(f[ix+3].split()[3])
#
#        # Now make the computer names uniform using predefined names
#        if self.host is not None:
#            for name in self._computer_names:
#                if name.lower() in self.host: self.host = name
