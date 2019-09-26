#!/usr/bin/env python

import numpy
from math import sqrt
import m8r

# initialize parameters
par = m8r.Par()

# input and output
vel=m8r.Input()  
tim=m8r.Output() 

# time axis from input
nt = vel.int('n1')
dt = vel.float('d1')

# offset axis from command line
nh = par.int('nh',1)      # number of offsets
dh = par.float('dh',0.01) # offset sampling
h0 = par.float('h0',0.0)  # first offset

# get reflectors
nr = par.int('nr',1) # number of reflectors
r = par.ints('r',nr)

type = par.string('type','hyperbolic')
# traveltime computation type

niter = par.int('niter',10)
# maximum number of shooting iterations

# put dimensions in output
tim.put('n1',nh)
tim.put('d1',dh)
tim.put('o1',h0)
tim.put('n2',nr)

# read velocity
v = numpy.zeros(nt,'f')
vel.read(v)

# convert to velocity squared
v = v*v

t = numpy.zeros(nh,'f')

for ir in range(nr):
    nt = r[ir]
    t0 = nt*dt # zero-offset time 
    t2 = t0*t0

    p = 0.0

    for ih in range(nh):
        h = h0+ih*dh # offset
        h2 = h*h

        if type[0] == 'h':
            # hyperbolic approximation 
            v2 = numpy.sum(v)/nt
            t[ih] = sqrt(t2+h2/v2)

        elif type[0] == 's':
            # shifted hyperbola 
         
            ### MODIFY BELOW ###

            s = 0.0
            v2 = numpy.sum(v)/nt
            t[ih] = sqrt(t2+h2/v2)

        elif type[0] == 'e': 
            # exact 
	    
            ### MODIFY BELOW ###

            for iter in range(niter): 
                hp = 0.0
                for it in range(nt):
                    v2 = v[it]
                    hp += v2/sqrt(1.0-p*p*v2)
                hp *= p*dt

            ### SOLVE h(p)=h ###

            tp = 0.0
            for it in range(nt):
                v2 = v[it]
                tp += dt/sqrt(1.0-p*p*v2)
            t[ih] = tp
        else:
            raise RuntimeError('Unknown type')

    tim.write(t)
