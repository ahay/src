#!/usr/bin/env python
'''
Multiple 1D inputs clipped to length of reference input

Built for log data manipulation - clips six input logs to the length of the reference well log
'''

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

par = rsf.Par()
depth_name = par.string('depth')
par = rsf.Par()
log1_name = par.string('log1')
par = rsf.Par()
log2_name = par.string('log2')
par = rsf.Par()
log3_name = par.string('log3')
par = rsf.Par()
log4_name = par.string('log4')

logrefa = rsf.Input() # Reference well log
deptha = rsf.Input(depth_name) # Depth information
log1a = rsf.Input(log1_name) # Well log 1
log2a = rsf.Input(log2_name) # Well log 2
log3a = rsf.Input(log3_name) # Well log 3
log4a = rsf.Input(log4_name) # Well log 4

assert 'float' == logrefa.type
assert 'float' == deptha.type
assert 'float' == log1a.type
assert 'float' == log2a.type
assert 'float' == log3a.type
assert 'float' == log4a.type

sdepth = deptha.int("n1")

logref = np.zeros(sdepth,'f')
depth = np.zeros(sdepth,'f')
log1 = np.zeros(sdepth,'f')
log2 = np.zeros(sdepth,'f')
log3 = np.zeros(sdepth,'f')
log4 = np.zeros(sdepth,'f')

logrefa.read(logref)
deptha.read(depth)
log1a.read(log1)
log2a.read(log2)
log3a.read(log3)
log4a.read(log4)

minsamp = 0
maxsamp = sdepth

for i in range(sdepth - 1):
    if ((logref[i] != 0) & (logref[i + 1] != 0)):
        break
    if ((logref[i] == 0) & (logref[i + 1] != 0)):
        minsamp = i + 1;
        break

for i in range(sdepth - 1):
    if ((logref[sdepth - i - 1] != 0) & (logref[sdepth - i - 2] != 0)):
        break
    if ((logref[sdepth - i - 1] == 0) & (logref[sdepth - i - 2] != 0)):
        maxsamp = sdepth - i - 2;
        break

logref_c = np.zeros(maxsamp - minsamp + 1)
depth_c = np.zeros(maxsamp - minsamp + 1)
log1_c = np.zeros(maxsamp - minsamp + 1)
log2_c = np.zeros(maxsamp - minsamp + 1)
log3_c = np.zeros(maxsamp - minsamp + 1)
log4_c = np.zeros(maxsamp - minsamp + 1)

for i in range(maxsamp - minsamp + 1):
    logref_c[i] = logref[minsamp + i]
    depth_c[i] = depth[minsamp + i]
    log1_c[i] = log1[minsamp + i]
    log2_c[i] = log2[minsamp + i]
    log3_c[i] = log3[minsamp + i]
    log4_c[i] = log4[minsamp + i]


logref_co = rsf.Output() # Output reference well log
depth_co = rsf.Output('depth_c') # Output depth information
log1_co = rsf.Output('log1_c') # Output clipped well log 1
log2_co = rsf.Output('log2_c') # Output clipped well log 2
log3_co = rsf.Output('log3_c') # Output clipped well log 3
log4_co = rsf.Output('log4_c') # Output clipped well log 4

logref_co.put('n1', maxsamp - minsamp + 1)
depth_co.put('n1', maxsamp - minsamp + 1)
log1_co.put('n1', maxsamp - minsamp + 1)
log2_co.put('n1', maxsamp - minsamp + 1)
log3_co.put('n1', maxsamp - minsamp + 1)
log4_co.put('n1', maxsamp - minsamp + 1)

logref_co.write(np.array(logref_c));
depth_co.write(np.array(depth_c));
log1_co.write(np.array(log1_c));
log2_co.write(np.array(log2_c));
log3_co.write(np.array(log3_c));
log4_co.write(np.array(log4_c));
