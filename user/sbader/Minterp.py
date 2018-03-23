#!/usr/bin/env python
'''
1D linear missing data interpolation

Linear interpolation of missing data (0 values) based on nearest nonzero samples
'''


import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

logrefa = rsf.Input() # Input 1D dataset

assert 'float' == logrefa.type
sdepth = logrefa.int("n1")

logref = np.zeros(sdepth,'f')
logrefa.read(logref)

min = 0
max = 0

out = np.zeros(sdepth)
for i in range(sdepth):
    out[i] = logref[i]

for i in range(sdepth):
    if (out[i] == 0):
        min = i - 1
        for j in range(sdepth - 1 - i):
            if (logref[i + j] != 0):
                max = i + j
                break
        for k in range(max - min):
            out[min + k] = logref[min] + (logref[max] - logref[min]) * k / (max - min)


logref_co = rsf.Output() # Output 1D dataset

logref_co.write(np.array(out));
