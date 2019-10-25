#!/usr/bin/env python
'''
Hard constraint 2D map

For use with LSIM if a specific alignment location is desired.
Use with reference and real datasets and scale LSIM scan by output
'''
from __future__ import division

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

# Inputs
par = rsf.Par()
logrefa = rsf.Input() # Input 1D data

assert 'float' == logrefa.type

sdepth = logrefa.int("n1")
sd = logrefa.float("d1")
so = logrefa.float("o1")

logref = np.zeros(sdepth,'f')
logrefa.read(logref)

num = par.float("value") # Location of hard constraint

window = par.int("wind") # Number of samples of hard constraint
assert window

out = np.zeros(sdepth)

loc = (num-so)/sd

if ( (loc < sdepth-window) &(loc > (window-1)) ):
    for i in range(window):
        out[int(loc-i)] = 1
        out[int(loc+i)] = 1

logref_co = rsf.Output()

logref_co.write(np.array(out)); # Output 1D data
