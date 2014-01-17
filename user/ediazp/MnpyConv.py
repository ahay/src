#!/usr/bin/env python
'''
Implements conv(a,b) along the fast axis 
a [file] : is taken from stdin
b [file] : is taken from  "flt"
Requires both files to have the same sampling interval

mode [string]:
 'full': returns an M+N-1 array, boundary effects are visible.
 'same': returns a max(M,N) array, boundary effects are visible.
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
mode = par.string("mode","same")

nc = {'full':na+nb-1,    'same':max(na,nb)}
oc = {'full':min(oa,ob), 'same':oa}

# output file
Fc = rsf.Output()
try:
  Fc.put("n1",nc[mode])
  Fc.put("o1",oc[mode])
except:
  print >> sys.stderr, \
      'mode %s'%mode+' is not recognized, valid modes are: "full" and "same"'
  sys.exit(1)

# ------------------------------------------------------------
n2 = Fa.size(1) # number of traces of input file
for i2 in range(n2):
  Fa.read(a)
  Fb.read(b)
  c = np.convolve(a,b,mode=mode)

  if norm:
    s = normalize(a,b)
    Fc.write(c*s)
  else:
    Fc.write(c)

# ------------------------------------------------------------
Fa.close()
Fb.close()
Fc.close()

