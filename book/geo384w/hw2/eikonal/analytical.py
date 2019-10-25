#!/usr/bin/env python

import numpy
from math import sqrt, hypot
import m8r

# initialize parameters
par = m8r.Par()

# input and output
inp=m8r.Input()  
out=m8r.Output()

# get grid dimensions
n1 = inp.int('n1')
n2 = inp.int('n2')
d1 = inp.float('d1')
d2 = inp.float('d2')

g1 = par.float('g1',0.0) # vertical gradient
g2 = par.float('g2',0.0) # horizontal gradient

gsq = g1*g1+g2*g2
g = sqrt(gsq)

v0 = par.float('v0')
# initial velocity or slowness squared

s = par.float('s',0.0)
# shot location at the surface

type = par.string('case','constant')
# case of velocity distribution

if 0.0 == g1 and 0.0 == g2:
    type='const'
    
time = numpy.zeros(n1,'f')

for i2 in range(n2):
   x2 = i2*d2
   for i1 in range(n1):
      x1 = i1*d1
      d = x1*x1+(x2-s)*(x2-s)
      
      if type[0] == 's':
      # slowness squared 
         s2 = v0+g1*x1+g2*x2
         z = 2.0*d/(s2+sqrt(s2*s2-gsq*d))
         time[i1] = (s2-gsq*z/6.0)*sqrt(z)

      elif type[0] == 'v':
      # velocity 
         s2 = 2.0*v0*(v0+g1*x1+g2*x2)

### CHANGE BELOW ### 
         time[i1] = hypot(x2-s,x1)/v0

      else:
         # constant velocity 
         time[i1] = hypot(x2-s,x1)/v0

   out.write(time)
