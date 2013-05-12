#!/usr/bin/env python
'''
This program implements Fornberg, 1988
paper for digital differentiators
of arbitrary order.

So, it computes first, second, n derivative.
'''


import rsf.api as rsf
import numpy as np


def Fornberg_filter(Nlenght,order):
  x = np.arange(-int(Nlenght/2),Nlenght/2+1,1)
  z = 0.0
  k = order
  n = Nlenght
  m = n -1
  c1 = 1
  c4 = x[0] - z
  c = np.zeros((n, m+1),'f')
  c[0,0] = 1
  for i in xrange(1, n):
    mn = min(i, m)
    c2 = 1
    c5 = c4
    c4 = x[i] - z
    for j in xrange(i):
      c3 = x[i] - x[j]
      c2 = c2*c3
      if j == i-1:
        for v in xrange(mn, 0, -1):
          c[i,v] = c1*(v*c[i-1,v-1] - c5*c[i-1,v])/c2
        c[i,0] = -c1*c5*c[i-1,0]/c2

      for v in xrange(mn, 0, -1):
        c[j,v] = (c4*c[j,v] - v*c[j,v-1])/c3
      c[j,0] = c4*c[j,0]/c3
    c1 = c2
  return c.T[k]  


par = rsf.Par()

Fa = rsf.Input()
Fout = rsf.Output()

# pars from command line
order  = par.int("order",1) # order of the derivative, default first derivative 
lenght = par.int("lenght",5) # filter lenght, the lenghtier the accurate, but also gets costlier 
scale  = par.bool("scale",True) # scales by 1/d^order
axis   = par.int("axis",1) # apply differentiator along axis, default is fast axis

# parameters from input file
n1 = Fa.int("n1")    # number of samples from input
d1 = Fa.float("d1")  # dt

a = np.zeros(n1,'f')
f = Fornberg_filter(lenght,order)

ntraces = Fa.size(1) # number of traces of input file
idn = 1./(d1)**order
for i2 in range(ntraces):
  Fa.read(a)
  out = np.convolve(f,a,'same')
  
  Fout.write(out)

Fa.close()
Fout.close()

