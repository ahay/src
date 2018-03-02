#!/usr/bin/env python
'''
make a checkerboard
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()
verb = par.bool('verb',False) # verbosity flag

N  = par.int( 'N',1)

nx = par.int  ('nx',1)
ox = par.float('ox',0.0)
dx = par.float('dx',1.0)

ny = par.int  ('ny',1)
oy = par.float('oy',0.0)
dy = par.float('dy',1.0)

nz = par.int  ('nz',1)
oz = par.float('oz',0.0)
dz = par.float('dz',1.0)

Fou = rsf.Output() # output file

Fou.put("n1",nz)
Fou.put("o1",oz)
Fou.put('d1',dz)

Fou.put("n2",nx)
Fou.put("o2",ox)
Fou.put("d2",dx)

Fou.put("n3",ny)
Fou.put("o3",oy)
Fou.put("d3",dy)

if(verb):
    print >>sys.stderr,N,nx,ny,nz
    
# ------------------------------------------------------------
ccc = np.zeros( (ny,nx,nz) ,'f')

cy = False
for iy in range(ny):
    if np.mod(iy,N)==0:
        cy = not cy;
        
    cx = True
    for ix in range(nx):
        if np.mod(ix,N)==0:
            cx = not cx;

        cz = True
        for iz in range(nz):
            if np.mod(iz,N)==0:
                cz = not cz;

            if( (    cy and (    cx and     cz)) or
                (    cy and (not cx and not cz)) or
                (not cy and (    cx and not cz)) or
                (not cy and (not cx and     cz)) ):
                ccc[iy][ix][iz]=1
# ------------------------------------------------------------

Fou.write(ccc)

Fou.close()
