#!/usr/bin/env python
'''
Residual well log realignment

Built for log data manipulation - algined axis (o1/d1/n1) with reference well log
'''

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

par = rsf.Par()
log1_name = par.string('log1')

logrefa = rsf.Input() # Reference Well Log
log1a = rsf.Input(log1_name) # Input Well Log

assert 'float' == logrefa.type
assert 'float' == log1a.type

sdepth = logrefa.int("n1")
ldepth = log1a.int("n1")
dlt = logrefa.float("d1")

olt = logrefa.float("o1")

logref = np.zeros(sdepth,'f')
log1 = np.zeros(ldepth,'f')

logrefa.read(logref)
log1a.read(log1)

minsamp = 0

maxs = sdepth - 1

if (logref[0] == 0):
    for i in range(maxs):
        if ((logref[i] != 0) & (logref[i + 1] != 0)):
            break
        if ((logref[i] == 0) & (logref[i + 1] != 0)):
            minsamp = i + 1;
            break


log1_co = rsf.Output()

if (minsamp != 0):
    log1_co.put('o1', minsamp*dlt + olt)


log1_co.put('n1', ldepth)
log1_co.write(np.array(log1));
