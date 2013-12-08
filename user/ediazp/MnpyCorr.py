#!/usr/bin/env python
'''
Implements corr(a,b) along the fast axis 
a [file] : is taken from stdin
b [file] : is taken from  "flt"
Requires both files to have the same sampling interval
'''

import rsf.api as rsf
import numpy as np
import sys

# ------------------------------------------------------------
def normalize(a,b):
  a0 = np.sqrt((a**2).sum())
  b0 = np.sqrt((b**2).sum())
  s = 1./(a0*b0)
  return s
# ------------------------------------------------------------

par = rsf.Par()

# File "a"
Fa = rsf.Input()
na = Fa.int  ("n1")
oa = Fa.float("o1")
da = Fa.float("d1")
a = np.zeros(na,'f')

# File "b"
Fb = rsf.Input("flt")
nb = Fb.int  ("n1")
ob = Fb.float("o1")
db = Fb.float("d1")
b = np.zeros(nb,'f')

if not da == db:
  print  >> sys.stderr, 'input files must have the same sampling interval'
  sys.exit(1)
  
# ------------------------------------------------------------

# command line params
norm = par.bool("norm",False) # normalize output
nc = par.int("nc",100)    # number of correlation lags 
if nc > na or nc > nb: nc = min(na-1,nb-1) 

# ------------------------------------------------------------
  
# output file
Fc = rsf.Output()
Fc.put("n1", nc*2+1)
Fc.put("o1",-nc*da)
Fc.put('d1',    da)

# ------------------------------------------------------------

# center sample of the correlation
center = (na+nb-1)/2 - (na-nb)*0.5
l1 = center-nc
l2 = center+nc+1

# ------------------------------------------------------------

n2 = Fa.size(1) # number of traces of input file
for i2 in range(n2):
  Fa.read(a)
  Fb.read(b)
  c = np.correlate(a,b,mode='full')
  
  if norm:
    s = normalize(a,b)
    Fc.write(c[l1:l2]*s)
  else:
    Fc.write(c[l1:l2])

# ------------------------------------------------------------
Fa.close()
Fb.close()
Fc.close()

