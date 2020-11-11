#! /usr/bin/env python

from __future__ import print_function

def progress_bar(current, last, bar='=', empty='-', sep='|'):
    '''Displays a progress bar on screen.
    NB: Nothing else can be printed between successive
    calls of the progress bar.'''

    import sys

    step = last / 100.
    pcnt = int(current/step)
    end = ''
    sp = sep
    if current == last:
        end = '\n'
        sp = ''
        pcnt = 100
    print ('{0}{1}{2}{3}{4}{5:3}%\r'.format(sep,  bar*pcnt,
           sp, empty*(99-pcnt), sep, pcnt), end=end)
    sys.stdout.flush()
