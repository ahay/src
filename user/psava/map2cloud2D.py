#!/usr/bin/env python
'''
compute normals of a 3D surface
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag

# ------------------------------------------------------------
Fin = rsf.Input()             # input file
nx  = Fin.int  ("n1")
ox  = Fin.float("o1")
dx  = Fin.float("d1")
#print  >> sys.stderr,nx,ox,dx

z   = np.zeros( nx,'f')
Fin.read(z)

# ------------------------------------------------------------
Fou = rsf.Output()            # output file
Fou.put("n1",4)
Fou.put("o1",0)
Fou.put('d1',1)

Fou.put("n2",nx)
Fou.put("o2",0)
Fou.put("d2",1)

dou = np.zeros(4,'f')

# ------------------------------------------------------------
for ix in range(nx):
    x = ox + ix * dx

    ax = dx
    if(ix < nx-1):
        az = z[ix+1] - z[ix]
    else:
        az = z[ix] - z[ix-1]

    cx = -az
    cz = +ax
    cc = np.sqrt(np.power(cx,2)+np.power(cz,2))

    dou = np.array([ x, z[ix], cx/cc,cz/cc])
    Fou.write(dou)

# ------------------------------------------------------------
Fin.close()
Fou.close()
