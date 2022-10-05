#!/usr/bin/env python

import sys
import math
import numpy as np
import m8r

# initialize
par = m8r.Par()
inp = m8r.Input()
out = m8r.Output()

# get dimensions from input
n1 = inp.int('n1')
n2 = inp.int('n2')

# get parameters from command line
angle = par.float('angle',90.)
# rotation angle

interp = par.string('interp','nearest')
# [n,l,c] interpolation type 

# convert degrees to radians
angle = angle*math.pi/180.
cosa = math.cos(angle)
sina = math.sin(angle)

orig = np.zeros([n1,n2],'f')
rotd = np.zeros([n1,n2],'f')

# read data
inp.read(orig)

# central point 
c1 = (n1-1)*0.5
c2 = (n2-1)*0.5

for i2 in range(n2):
    for i1 in range(n1):
        # rotated coordinates
        x1 = c1+(i1-c1)*cosa-(i2-c2)*sina
        x2 = c2+(i1-c1)*sina+(i2-c2)*cosa

	    # nearest neighbor
        k1 = int(math.floor(x1))
        k2 = int(math.floor(x2))
        x1 -= k1	    
        x2 -= k2

        if interp[0] == 'n':
            # nearest neighbor
            if x1 > 0.5:
                k1 += 1
            if x2 > 0.5:
                k2 += 1
            if k1 >=0 and k1 < n1 \
              and k2 >=0 and k2 < n2:
                rotd[i2,i1] = orig[k2,k1]
            else:
                rotd[i2,i1] = 0.
        elif interp[0] == 'l':
            # bilinear
            if k1 >=0 and k1 < n1-1 \
              and k2 >=0 and k2 < n2-1:
                rotd[i2,i1] = \
			    (1.-x1)*(1.-x2)*orig[k2,k1]   + \
			    x1     *(1.-x2)*orig[k2,k1+1] + \
			    (1.-x1)*x2     *orig[k2+1,k1] + \
			    x1     *x2     *orig[k2+1,k1+1]
            else:
                rotd[i2,i1] = 0.
        elif interp[0] == 'c':
            # cubic convolution 
		    # !!! ADD CODE !!!
            break
        else:
            sys.stderr.write('Unknown method "%s"' % interp)
            sys.exit(1)

 # write result */
out.write(rotd)
sys.exit(0)
