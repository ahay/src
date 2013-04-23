#!/usr/bin/env python
'''
This program is a wrapper for 
numpy.correlate . It uses the 'full' mode
correlation. 

It implements corr(a,b) along the fast axis 

a [file] : is taken from stdin
b [file] : is taken from  signal

As for now, it requires both files to have the 
same sampling interval
'''


import rsf.api as rsf
import numpy as np
import sys

def normalize(a,b):
  a0 = np.sqrt((a**2).sum())
  b0 = np.sqrt((b**2).sum())
  s = 1./(a0*b0)

  return s



par = rsf.Par()

Fa = rsf.Input()
Fb = rsf.Input("signal")
Fout = rsf.Output()

# pars from command line
ntau = par.int("ntau",100) # number of tau lags 
norm = par.bool("norm",False) # if True scale correlation by 1/(sqrt(corr(a,a)[0])*sqrt(corr(b,b)[0])) so the correlation goes from -1,1
print norm



# parameters from input file
n1a = Fa.int("n1")    # number of samples from input
n1b = Fb.int("n1")    # number of samples from input
d1a = Fa.float("d1")  # dt
d1b = Fb.float("d1")  # dt

if not d1a == d1b:
  print  >> sys.stderr, 'both files are required to have the same sampling interval'
  sys.exit(1)


if ntau > n1a or ntau > n1b: ntau = min(n1a-1,n1b-1) 

ntraces = Fa.size(1) # number of traces of input file

a = np.zeros(n1a,'f')
b = np.zeros(n1b,'f')
center = (n1a+n1b-1)/2 -(n1a-n1b)*0.5 # center sample of the correlation
l1 = center-ntau
l2 = center+ntau+1

Fout.put("n1",ntau*2+1)
Fout.put("o1",-ntau*d1a)

for i2 in range(ntraces):
  Fa.read(a)
  Fb.read(b)
  corr = np.correlate(a,b,mode='full')
  out = corr[l1:l2]
  if norm:
    s = normalize(a,b)
    Fout.write(out*s)
  else:
    Fout.write(out)

Fa.close()
Fb.close()
Fout.close()

