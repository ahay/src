#!/usr/bin/env python
'''
output point cloud from gridded surface
-> cloud components: x,y,z, nx,ny,nz, vx,vy,vz
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
sphc = par.bool('sphc',False) # spherical coordinates flag
# print(sphc, file=sys.stderr)
deg2rad = np.pi / 180         # degrees to radians

# ------------------------------------------------------------
Fin = rsf.Input()             # input file
n1  = Fin.int  ("n1")
o1  = Fin.float("o1")
d1  = Fin.float("d1")

n2  = Fin.int  ("n2")
o2  = Fin.float("o2")
d2  = Fin.float("d2")

r   = np.zeros( (n2,n1),'f')  # read elevation
Fin.read(r)

# ------------------------------------------------------------
Fou = rsf.Output()            # output file
Fou.put("n1",9)
Fou.put("o1",0)
Fou.put('d1',1)

Fou.put("n2",n2*n1)
Fou.put("o2",0)
Fou.put('d2',1)

dou = np.zeros(9,'f')

# ------------------------------------------------------------
# compute Cartesian coordinates
x = np.zeros( (n2,n1),'f')
y = np.zeros( (n2,n1),'f')
z = np.zeros( (n2,n1),'f')

if sphc: # spherical to Cartesian coordinates
    for i2 in range(n2):
        lon = o2 + i2 * d2 # [deg] longitude
        lon *= deg2rad     # [rad] longitude

        for i1 in range(n1):
            lat = o1 + i1 * d1 # [deg] latitude
            lat *= deg2rad     # [rad] latitude

            x[i2,i1] = np.cos(lat) * np.cos(lon)
            y[i2,i1] = np.cos(lat) * np.sin(lon)
            z[i2,i1] = np.sin(lat)

    # scale by elevation
    x *= r
    y *= r
    z *= r

else: # local Cartesian coordinates
    for i2 in range(n2):
        for i1 in range(n1):
            x[i2,i1] = o2 + i2 * d2
            y[i2,i1] = o1 + i1 * d1
            z[i2,i1] = r[i2,i1]

# ------------------------------------------------------------
# compute normals
for i2 in range(n2):
    for i1 in range(n1):

        # vector a
        if(i2 < n2-1):
            ax = x[i2+1][i1] - x[i2  ][i1]
            ay = y[i2+1][i1] - y[i2  ][i1]
            az = z[i2+1][i1] - z[i2  ][i1]
        else:
            ax = x[i2  ][i1] - x[i2-1][i1]
            ay = y[i2  ][i1] - y[i2-1][i1]
            az = z[i2  ][i1] - z[i2-1][i1]

        # vector b
        if(i1 < n1-1):
            bx = x[i2][i1+1] - x[i2][i1  ]
            by = y[i2][i1+1] - y[i2][i1  ]
            bz = z[i2][i1+1] - z[i2][i1  ]
        else:
            bx = x[i2][i1  ] - x[i2][i1-1]
            by = y[i2][i1  ] - y[i2][i1-1]
            bz = z[i2][i1  ] - z[i2][i1-1]

        # normal vector n = a x b
        nx = ay*bz - by*az
        ny = az*bx - bz*ax
        nz = ax*by - bx*ay
        nn = np.sqrt(np.power(nx,2) +
                     np.power(ny,2) +
                     np.power(nz,2) )

        # output x,y,z, nx,ny,nz
        dou = np.array([ x[i2][i1],  y[i2][i1],  z[i2][i1],
                       -nx/nn,     -ny/nn,     -nz/nn,
                         0,          0,          0 ])
        Fou.write(dou)

# ------------------------------------------------------------
Fin.close()
Fou.close()
