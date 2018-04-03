#!/usr/bin/env python
'''
Estimate energy of input

E(t) = \sum\limits_{k=(t-\frac{R}{2})}^{(t+\frac{R}{2})}A(k)^2
'''


import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

# Inputs
par = rsf.Par()
input = rsf.Input() # Input 1D dataset

assert 'float' == input.type
sdepth = input.int("n1")

log = np.zeros(sdepth,'f')
input.read(log)

log_eo = rsf.Output() # Output 1D energy calculation

radius = par.int("wind") # Rolling window size
assert radius

log_e = np.empty(sdepth)

beg = 0
end = 0
maxv = 0

for i in range(radius):
    beg = beg + log[i]*log[i]

for i in range(radius):
    end = end + log[sdepth - i - 1]*log[sdepth - i - 1]

for i in range(sdepth):
    if (i < radius):
        log_e[i] = beg
    if (i > (sdepth - radius - 1)):
        log_e[i] = end
    else:
        middle = 0
        for j in range(2*radius):
            middle = middle + log[i - radius + j]*log[i - radius + j]
        log_e[i] = middle
            
    if (i > 0):
        if (log_e[i] > maxv):
            maxv = log_e[i]

for i in range(sdepth):
    log_e[i] = log_e[i]/maxv

log_eo.write(np.array(log_e));
