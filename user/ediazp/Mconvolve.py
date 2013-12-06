#!/usr/bin/env python
'''
This program is a wrapper for 
numpy.convolve . It uses the 'full' mode
correlation. 

It implements conv(a,b) along the fast axis 

a [file] : is taken from stdin
b [file] : is taken from  filter

mode [string]:
  'full': returns an M+N-1 array, boundary effects 
          are visible.

  'same': returns a max(M,N) array, boundary effects
          are visible.

As for now, it requires both files to have the 
same sampling interval
'''


import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

Fa = rsf.Input()
Fb = rsf.Input("filter")
Fout = rsf.Output()

# pars from command line
mode = par.string("mode","same")
print mode




# parameters from input file
n1a = Fa.int("n1")    # number of samples from input
n1b = Fb.int("n1")    # number of samples from input
d1a = Fa.float("d1")  # dt
d1b = Fb.float("d1")  # dt
o1a = Fa.float("o1")  # dt
o1b = Fb.float("o1") 

length = {'full':n1a+n1b-1,
          'same':max(n1a,n1b)}

if not d1a == d1b:
  print  >> sys.stderr, \
           'both files are required to have the same sampling interval'
  sys.exit(1)


ntraces = Fa.size(1) # number of traces of input file

a = np.zeros(n1a,'f')
b = np.zeros(n1b,'f')


origin = {'full':min(o1a,o1b),
          'same':o1a}
try:
  Fout.put("n1",length[mode])
  Fout.put("o1",origin[mode])
except:
  print >> sys.stderr, \
      'mode %s'%mode+' is not recognized, valid modes are: "full" and "valid"'
  sys.exit(1)

for i2 in range(ntraces):
  Fa.read(a)
  Fb.read(b)
  conv = np.convolve(a,b,mode=mode)
  Fout.write(conv)

Fa.close()
Fb.close()
Fout.close()

