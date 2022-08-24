#!/usr/bin/env python

import sys
import math
import numpy
import m8r

# initialization
par = m8r.Par()
inp = m8r.Input()
out = m8r.Output()

# get trace parameters
nt = inp.int('n1')
dt = inp.float('d1')
t0 = inp.float('o1')
 
#  get number of traces
nx = inp.leftsize(1)

na = par.int('na',1)     # number of alpha values
da = par.float('da',0.0) # increment in alpha
a0 = par.float('a0',0.0) # first value of alpha

# change output data dimensions
out.put('n1',na)
out.put('n2',1)
out.put('d1',da)
out.put('o1',a0)

trace = numpy.zeros(nt,'f')
tgain  = numpy.zeros(nt,'f')
ofunc = numpy.zeros(na,'f')

# loop over traces
for ix in range(nx):
    # read data
    inp.read(trace)

    # loop over alpha
    for ia in range(na):
        alpha = a0+ia*da

        # loop over time samples
        for it in range(nt):
            t = t0+it*dt

            # apply gain t^alpha 
            s = trace[it]*math.pow(t,alpha)
		
            # !!! MODIFY THE NEXT LINE !!! 
            ofunc[ia] += s*s

# write output
out.write(ofunc)
sys.exit(0)
