#!/usr/bin/env python
'''
Inverse Multiquadratic Radial Basis Function

1/sqrt(1 + (eps*r)^2) where r is distance from source
'''


import sys
import numpy as np
import rsf.api as rsf
import string
import math as mt

# Inputs
par = rsf.Par()
cubea = rsf.Input() # Seismic Cube (Used for size of output)

assert 'float' == cubea.type

t_dim = cubea.int("n1")
xl_dim = cubea.int("n2")
il_dim = cubea.int("n3")

rbf_out = rsf.Output() # Output RBF cube

cube = cubea.read()

xl_pos = par.int("xl") # n2 location of source
il_pos = par.int("il") # n3 location of source
eps = par.float("eps") # Scalar factor

bound = par.int("boundary") # Scalar factor
if (bound == 1):
    otherdata = par.string('other') # Boundary map
    othera = rsf.Input(otherdata)
    assert 'float' == othera.type

    other = othera.read()


rbf=np.zeros((il_dim,xl_dim, t_dim))

if (bound == 0):
    for i in range(il_dim):
        for j in range(xl_dim):
            rbf[i][j][:] = 1/mt.pow(1+mt.pow(eps*mt.pow(mt.pow(i-il_pos,2)+mt.pow(j-xl_pos,2),0.5),2), 0.5)
else:
    for i in range(il_dim):
        for j in range(xl_dim):
            rbf[i][j][:] = other[i][j]/mt.pow(1+mt.pow(eps*mt.pow(mt.pow(i-il_pos,2)+mt.pow(j-xl_pos,2),0.5),2), 0.5)

rbf_out.write(np.array(rbf));

rbf_out.close()
