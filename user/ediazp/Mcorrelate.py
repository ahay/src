#!/usr/bin/env python
'''
This program is a wrapper for 
numpy.correlate . It uses the 'full' mode
correlation. 

It implements corr(a,b) along the fast axis 

a [file] : is taken from stdin
b [file] : is taken from  signal

'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

Fa = rsf.Input()
Fb = rsf.Input("signal")
Fout = rsf.Output()

# pars from command line
ntau = par.int("ntau",100) # number of tau lags 

# parameters from input file
n1a = Fa.int("n1")    # number of samples from input
n1b = Fb.int("n1")    # number of samples from signal

d1 = Fa.float("d1")  # dt

if ntau > n1a or ntau > n1b : ntau = min(n1a-1,n1b-1)

ntraces = Fa.size(1) # number of traces of input file



#######################################################################
# Check that files have the same sampling interval and geometry
if not ntraces == Fb.size(1):
  sys.stderr.write('error! :input and signal do not have the same number\
                    of traces\n')
  sys.exit(1)

if not d1 == Fb.float("d1"):
  sys.stderr.write('error! :input and signal do not have the same \
                    sampling interval\n')
  sys.exit(1)
#######################################################################
  

  

a = np.zeros(n1a,'f')
b = np.zeros(n1b,'f')

center = (n1a+n1b-1)/2 # center sample of the correlation
l1 = center-ntau
l2 = center+ntau+1

Fout.put("n1",ntau*2+1)
Fout.put("o1",-ntau*d1)


for i2 in range(ntraces):
  Fa.read(a)
  Fb.read(b)
  corr = np.correlate(a,b,mode='full')
  out = corr[l1:l2]
  Fout.write(out)

Fa.close()
Fb.close()
Fout.close()

