#!/usr/bin/env python
'''
One-sided data clipping

Built for log data manipulation - remove extraneous values introduced from LSIM shifting
'''


import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

# Inputs
par = rsf.Par()

logrefa = rsf.Input() # Input 1D trace or log

assert 'float' == logrefa.type
sdepth = logrefa.int("n1")

logref = np.zeros(sdepth,'f')
logrefa.read(logref)

num = par.float("value") # Output if log is less than clip

clip = par.float("clip") # Clipping value
assert clip

out = np.zeros(sdepth)

for i in range(sdepth):
    if (logref[i] > clip):
        out[i] = logref[i]
    else:
        out[i] = num

logref_co = rsf.Output() # Output 1D trace or log

logref_co.write(np.array(out));
