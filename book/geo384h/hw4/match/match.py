#!/usr/bin/env python

import sys
import numpy as np
import m8r

# initialize
par = m8r.Par()
inp = m8r.Input()
out = m8r.Output()
oth = m8r.Input('other')

adj = par.bool('adj',False) # adjoint flag 

if adj:
    # input data, output filter
    n1 = inp.int('n1')
    n2 = inp.int('n2')
    nf = par.int('nf') # filter size

    out.put('n1',nf)
    out.put('n2',1)
else:
	# input filter, output data
    nf = inp.int('n1')
    n1 = oth.int('n1')
    n2 = oth.int('n2')

    out.put('n1',n1)
    out.put('n2',n2)
    
filt = np.zeros(nf,'f')
data = np.zeros(n1,'f')
noiz = np.zeros(n1,'f')

if not adj:
    inp.read(filt)

for i2 in range(n2):
    oth.read(noiz)

    if adj:
	    inp.read # !!! COMPLETE LINE !!
    else:
        data[:] = 0.

    for i in range(nf):
        for i1 in range(n1):
            j=i1-i+nf//2 # symmetric filter 

            # zero value boundary conditions 
            if j < 0 or j >= n1:
                continue
		    
            if adj:
                filt[i] += # !!! COMPLETE LINE !!! 
            else:
                data[i1] += noiz[j]*filt[i]

    if not adj:
        out.write(data)

if adj:
    out.write  # !!! COMPLETE LINE !!!
		 
sys.exit(0)
