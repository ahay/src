#!/usr/bin/env python
'''
Dataset padding/clipping

Built for log data manipulation - Before and after LSIM alignments to deal with edge-effects and extraneous values introduced from LSIM shifting
'''


import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

par = rsf.Par()
loga = rsf.Input() # Input well log

assert 'float' == loga.type
sdepth = loga.int("n1")
olog = loga.float("o1")

log = np.zeros(sdepth,'f')
loga.read(log)

log_eo = rsf.Output() # Output well log

switch = par.int("switch") # (0 = Two-sided axis extension by first and last non-zero sample in dataset); (2 = Two-sided axis reduction); (3 = Matches starting value and number of samples between input and reference well log); (else = pad data to dataset size by first and last nonzero sample); (4 = Testing)


if ( (switch == 0) | (switch == 2) | (switch == 4)):
    num = par.int("val") # Sample manipulation (switch=0/2)
    assert num

if ( switch == 3 ):
    reflog = par.string('reflog') # Reference log (switch=3)
    refa = rsf.Input(reflog)
    assert 'float' == refa.type
    oref = refa.float("o1")
    nref = refa.int("n1")

    ref = np.zeros(nref,'f')
    refa.read(ref)

minsamp = 0
maxsamp = sdepth - 1

if (log[0] == 0):
    for i in range(sdepth - 1):
        if ((log[i] != 0) & (log[i + 1] != 0)):
            break
        if ((log[i] == 0) & (log[i + 1] != 0)):
            minsamp = i + 1;
            break

if (log[sdepth - 1] == 0):
    for i in range(sdepth - 1):
        if ((log[sdepth - i - 1] != 0) & (log[sdepth - i - 2] != 0)):
            break
        if ((log[sdepth - i - 1] == 0) & (log[sdepth - i - 2] != 0)):
            maxsamp = sdepth - i - 2;
            break

if (switch == 0):
    minval = log[minsamp + 50]
    maxval = log[maxsamp - 50]

    log_a = np.empty(minsamp + 50 + num)
    log_a.fill(minval)

    log_c = np.empty(sdepth - maxsamp + 50 + num)
    log_c.fill(maxval)

    log_e = np.concatenate((log_a, log[(minsamp + 50):(maxsamp - 50)], log_c), axis=0)
    log_eo.put('n1', sdepth + 2*num)

elif (switch == 2):
    log_e = log[num:(sdepth-num)]
    log_eo.put('n1', sdepth - 2*num)

elif (switch == 3):
    minval = log[minsamp]
    maxval = log[maxsamp]

    if (oref < olog):

        log_a = np.empty(2*(olog-oref) + minsamp)
        log_a.fill(minval)
    
        log_c = np.empty(nref - 2*(olog-oref) - (maxsamp - minsamp) )
        log_c.fill(maxval)

        log_e = np.concatenate((log_a, log[minsamp:maxsamp], log_c), axis=0)
        log_eo.put('o1', oref)
        log_eo.put('n1', nref)

    elif ( (oref == olog) & (sdepth < nref) ):
        log_c = np.empty(nref - (maxsamp - minsamp) )
        log_c.fill(maxval)

        log_e = np.concatenate((log, log_c), axis=0)
        log_eo.put('n1', nref)

    else:
        log_e = log

elif (switch == 4):
    minval = log[0]
    maxval = log[sdepth - 1]

    log_a = np.empty(num)
    log_a.fill(minval)

    log_c = np.empty(num)
    log_c.fill(maxval)

    log_e = np.concatenate((log_a, log, log_c), axis=0)
    log_eo.put('n1', sdepth + 2*num)

else:
    minval = log[minsamp]
    maxval = log[maxsamp]

    log_a = np.empty(minsamp)
    log_a.fill(minval)

    log_c = np.empty(sdepth - maxsamp)
    log_c.fill(maxval)

    log_e = np.concatenate((log_a, log[(minsamp):(maxsamp)], log_c), axis=0)




log_eo.write(np.array(log_e));
