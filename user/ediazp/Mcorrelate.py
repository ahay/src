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


par = rsf.Par()

b = par.string("signal",None) # File name with other signal

ntau = par.int("ntau",100) # number of tau lags 


Fa = rsf.Input()
Fb = rsf.Input(b)
Fout = rsf.Output()



n1 = Fa.int("n1")    # number of samples
d1 = Fa.float("d1")  # dt

if ntau > n1 : ntau=n1-1

ntraces = Fa.size(1) # number of traces

a = np.zeros(n1,'f')
b = np.zeros(n1,'f')

center = (n1*2-1)/2 # center sample of the correlation
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

