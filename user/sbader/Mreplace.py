#!/usr/bin/env python
'''
1D dataset padding 

Built for time-depth relationship manipulation. LSIM scan/pick is smooth intoducing non-real updates at top of pick.
'''

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

par = rsf.Par()
loga = rsf.Input() # Input 1D dataset

assert 'float' == loga.type
sdepth = loga.int("n1")

log = np.zeros(sdepth,'f')
loga.read(log)

loc = par.int("loc") # Location of value used to replace
assert loc
num = par.int("num") # Number of values to replace at beginning of dataset
assert num

log_eo = rsf.Output() # Output 1D dataset

log_e = np.empty(sdepth)

for i in range(sdepth):
    if (i < num):
        log_e[i] = log[loc]
    else:
        log_e[i] = log[i]

log_eo.write(np.array(log_e));
