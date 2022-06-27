#!/usr/bin/env python
'''
time shift in the frequency domain
'''
import rsf.api as rsf
import numpy as np

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag

# ------------------------------------------------------------
Fin = rsf.Input()             # input file
nx = Fin.int  ("n1")
ox = Fin.float("o1")
dx = Fin.float("d1")

nf = Fin.int  ("n2")
of = Fin.float("o2")
df = Fin.float("d2")
Fou = rsf.Output()            # output file

# ------------------------------------------------------------
Ftt = rsf.Input('tt')
tt = np.zeros(nx,'f')
Ftt.read(tt)
Ftt.close()
# ------------------------------------------------------------

din = np.zeros(nx,dtype='complex64')
dou = np.zeros(nx,dtype='complex64')

for jf in range(nf):
    f = of + jf * df

    Fin.read(din)

    dou = din * np.exp(2j * np.pi * f * tt)

    Fou.write(dou)

# ------------------------------------------------------------
Fin.close()
Fou.close()
