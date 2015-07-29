#!/usr/bin/env python
'''
This program implements Fornberg, 1988
paper for digital differentiators
of arbitrary order.

So, it computes first, second, n derivative along axis 1,2 or 3.
'''


import rsf.api as rsf
import numpy as np


def Fornberg_filter(Nlength,order):
  '''
  This function computes the Fornberg filters
  '''

  x = np.arange(-int(Nlength/2),Nlength/2+1,1)
  z = 0.0
  k = order
  n = Nlength
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


def goFilter1(Fin,Fout,f,order,scale=True):
  n1 = Fin.int("n1")
  d1 = Fin.float("d1")
  ntraces = Fin.size(1)

  idn = np.ones(1,'f')
  idn[0] = (1./(d1**order))

  a = np.zeros(n1,'f')

  for ir in range(ntraces):
    Fin.read(a)
    out = np.convolve(f,a,'same')
    if scale:
      out *= idn
    Fout.write(out)

def goFilter2(Fin,Fout,f,order,scale=True):
  n1 = Fin.int("n1")
  n2 = Fin.int("n2")
  d2 = Fin.float("d2")
  
  ntraces = Fin.size(2)
  idn = np.ones(1,'f')
  idn[0] = (1./(d2**order))
  a = np.zeros((n2,n1),'f')
  aux = np.zeros(n2,'f')

  for i3 in range(ntraces):
    Fin.read(a)
    for i1 in range(n1):
      aux = a[:,i1]
      out = np.convolve(f,aux,'same')
      a[:,i1] = out
    if scale:
      a *= idn
    Fout.write(a)
    
def goFilter3(Fin,Fout,f,order,scale=True):
  n1 = Fin.int("n1")
  n2 = Fin.int("n2")
  n3 = Fin.int("n3")
  d3 = Fin.float("d3")
  
  ntraces = Fin.size(3)
  idn = np.ones(1,'f')
  idn[0] = (1./(d3**order))

  a = np.zeros((n3,n2,n1),'f')
  aux = np.zeros(n3,'f')

  for i4 in range(ntraces):
    Fin.read(a)
    for i2 in range(n2):
      for i1 in range(n1):
        aux = a[:,i2,i1]
        out = np.convolve(f,aux,'same')
        a[:,i2,i1] = out
    if scale:
      a *= idn
    Fout.write(a)



#############################################################################################
#                               MAIN PROGRAM BEGINS:                                        #
#############################################################################################
par = rsf.Par()
Fin = rsf.Input()
Fout = rsf.Output()

functions = {1:goFilter1, 2:goFilter2, 3:goFilter3}
# pars from command line
order  = par.int("order",1) # order of the derivative, default first derivative 
length = par.int("length",5) # filter length, the lengthier the accurate, but also gets costlier 
scale  = par.bool("scale",True) # scales by 1/d^order
axis   = par.int("axis",1) # apply differentiator along axis, default is fast axis

f = Fornberg_filter(length,order)
try:
  functions[axis](Fin,Fout,f,order,scale)
except:
  import sys
  print >> sys.stderr, '========== sfnderiv ERROR ==========='
  print >> sys.stderr, 'Error: valid axis values are 1,2 or 3'
  print >> sys.stderr, '====================================='
  sys.exit(1)

Fin.close()
Fout.close()
