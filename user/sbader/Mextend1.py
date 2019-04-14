#!/usr/bin/env python
'''
Dataset padding - maintains dataset dims.

Pads dataset by first and last nonzero sample. This helps reduce artifacts introduced by PP well log interpolation.
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

log = np.zeros(sdepth,'f')
loga.read(log)

log_eo = rsf.Output() # Output well log

num = par.int("num") # Number of samples
assert num

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

minval = log[minsamp]
maxval = log[maxsamp]

maxnum = maxsamp
if ( (maxsamp+num) > (sdepth - 1) ):
    maxnum = sdepth - 1 - num

minnum = minsamp
if ( (minsamp-num) < 0 ):
    minnum = num

log_e = np.zeros(sdepth)
log_e[(minnum-num):minsamp] = minval
log_e[maxsamp:(maxnum+num)] = maxval
log_e[minsamp:maxsamp] = log[minsamp:maxsamp]

log_eo.write(np.array(log_e));
