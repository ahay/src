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
ny  = Fin.int  ("n1")
oy  = Fin.float("o1")
dy  = Fin.float("d1")
#print  >> sys.stderr,ny,oy,dy

nx  = Fin.int  ("n2")
ox  = Fin.float("o2")
dx  = Fin.float("d2")
#print  >> sys.stderr,nx,ox,dx

z   = np.zeros( (nx,ny),'f')
Fin.read(z)

# ------------------------------------------------------------
Fou = rsf.Output()            # output file
Fou.put("n1",6)
Fou.put("o1",0)
Fou.put('d1',1)

Fou.put("n2",nx*ny)
Fou.put("o2",0)
Fou.put('d2',1)

dou = np.zeros(6,'f')

# ------------------------------------------------------------
for ix in range(nx):
    x = ox + ix * dx

    for iy in range(ny):
        y = oy + iy * dy

        ax = dx
        ay = 0.0
        if(ix < nx-1):
            az = z[ix+1][iy] - z[ix][iy]
        else:
            az = z[ix][iy] - z[ix-1][iy]

        bx = 0.0
        by = dy
        if(iy < ny-1):
            bz = z[ix][iy+1] - z[ix][iy]
        else:
            bz = z[ix][iy] - z[ix][iy-1]

        cx = ay*bz - by*az
        cy = az*bx - bz*ax
        cz = ax*by - bx*ay
        cc = np.sqrt(np.power(cx,2)+np.power(cy,2)+np.power(cz,2))

        dou = np.array([ x, y, z[ix][iy], cx/cc,cy/cc,cz/cc])
        Fou.write(dou)

# ------------------------------------------------------------
Fin.close()
Fou.close()
