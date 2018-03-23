#!/usr/bin/env python
'''
1D dataset windowing

Built for log data manipulation - Clips padded values at the beginning and end of well log data to zero
'''

import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

par = rsf.Par()
loga = rsf.Input() # Input well log

log_eo = rsf.Output() # Output well log

switch = 0
switch = par.int("switch")

assert 'float' == loga.type
sdepth = loga.int("n1")

log = np.zeros(sdepth,'f')
loga.read(log)

minsamp = 0
maxsamp = sdepth - 1

avg_np = np.array(log)
avg = avg_np[avg_np > 0].mean()
if (avg < 10):
    rnd = 3
else:
    rnd = 0

minexit = 0
maxexit = 0

for i in range(15, sdepth - 16):
    if (minexit == 0):
        a = round(log[i], rnd)
        b = round(log[i+5], rnd)
        c = round(log[i+10], rnd)
        d = round(log[i+15], rnd)

    if ( (a > 0) & (a != b) & (b != c) & (c != d)):
        minsamp = i
        a = 0
        b = 0
        c = 0
        d = 0
        minexit = 1

    if (maxexit == 0):
        w = round(log[sdepth - i], rnd)
        x = round(log[sdepth - (i+5)], rnd)
        y = round(log[sdepth - (i+10)], rnd)
        z = round(log[sdepth - (i+15)], rnd)

    if ( (w > 0) & (w != x) & (x != y) & (y != z)):
        maxsamp = sdepth - i
        w = 0
        x = 0
        y = 0
        z = 0
        maxexit = 1

    if ( (minexit == 1) & (maxexit == 1) ):
        break

log_a = np.empty(minsamp)
log_c = np.empty(sdepth - maxsamp)

log_a.fill(0)
log_c.fill(0)

log_e = np.concatenate((log_a, log[(minsamp):(maxsamp)], log_c), axis=0)

log_eo.write(np.array(log_e))
