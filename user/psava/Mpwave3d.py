#!/usr/bin/env python
'''
generate a plane in t-x-y-z space
cx x + cy y + cz z - vel t = 0
(the plane goes through the origin)
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()
verb = par.bool('verb',False) # verbosity flag

nt=par.int("nt"); ot=par.float("ot"); dt=par.float("dt")
nx=par.int("nx"); ox=par.float("ox"); dx=par.float("dx")
ny=par.int("ny"); oy=par.float("oy"); dy=par.float("dy")
nz=par.int("nz"); oz=par.float("oz"); dz=par.float("dz")

vel= par.float("vel",1.0)
cx = par.float("cx",1.0)
cy = par.float("cy",1.0)
cz = par.float("cz",1.0)

if verb:
    print  >> sys.stderr,'t',nt,ot,dt,vel
    print  >> sys.stderr,'x',nx,ox,dx,cx
    print  >> sys.stderr,'y',ny,oy,dy,cy
    print  >> sys.stderr,'z',nz,oz,dz,cz

# ------------------------------------------------------------
Fdat = rsf.Output() 
Fdat.put("n1",nt); Fdat.put("o1",ot); Fdat.put('d1',dt)
Fdat.put("n2",nx); Fdat.put("o2",ox); Fdat.put('d2',dx)
Fdat.put("n3",ny); Fdat.put("o3",oy); Fdat.put('d3',dy)    
Fdat.put("n4",nz); Fdat.put("o4",oz); Fdat.put('d4',dz)    

dat = np.zeros(nt,'f')

# ------------------------------------------------------------
for iz in range(nz): 
    z = oz + iz*dz

    for iy in range(ny): 
        y = oy + iy*dy

        for ix in range(nx):
            x = ox + ix*dx

            for it in range(nt):
                dat[it]=0;        

            t = (cx/vel)*x + (cz/vel)*z
            it=int((t-ot)/dt)

            if it>0 and it<nt:
                dat[it]=1;
        
            Fdat.write(dat)

# ------------------------------------------------------------
Fdat.close()

