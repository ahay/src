#!/usr/bin/env python
'''
Oriented correlation
wfl [file] : is taken from stdin
opr [file] : is taken from  "opr"
Requires both files to have the same dimensions
'''
import rsf.api as rsf
import numpy as np
import sys

par = rsf.Par()

verb = par.bool('verb',False) # verbosity flag
adj  = par.bool( 'adj',False) # adjoint flag

Fwfl = rsf.Input()
Fopr = rsf.Input("opr")

nt=Fopr.int("n1"); ot=Fopr.float("o1"); dt=Fopr.float("d1")
nx=Fopr.int("n2"); ox=Fopr.float("o2"); dx=Fopr.float("d2")
ny=Fopr.int("n3"); oy=Fopr.float("o3"); dy=Fopr.float("d3")
nz=Fopr.int("n4"); oz=Fopr.float("o4"); dz=Fopr.float("d4")

# center of oriented correlation
ocox = par.float("ocox",0.0)
ocoy = par.float("ocoy",0.0)
ocoz = par.float("ocoz",0.0)

ixmid=(ocox-ox)/dx; 
iymid=(ocoy-oy)/dy; 
izmid=(ocoz-oz)/dz; 

if(nx==1): ixmid=0;
if(ny==1): iymid=0;
if(nz==1): izmid=0;

# ------------------------------------------------------------

opr = np.zeros((nz,ny,nx,nt),'f') # allocate opr array
wfl = np.zeros((nz,ny,nx,nt),'f') # allocate wfl array

Fopr.read(opr)              #     read opr array
Fwfl.read(wfl)              #     read wfl array

# ------------------------------------------------------------

ntlag = par.int("ntlag",100)
nxlag = par.int("nxlag",0)
nylag = par.int("nylag",0)
nzlag = par.int("nzlag",0)

nct=2*ntlag+1; oct=-ntlag*dt; dct=dt
ncx=2*nxlag+1; ocx=-nxlag*dx; dcx=dx
ncy=2*nylag+1; ocy=-nylag*dy; dcy=dy
ncz=2*nzlag+1; ocz=-nzlag*dx; dcz=dz

# ------------------------------------------------------------
Fcor = rsf.Output()        # output oriented cor

Fcor.put("n1",nct); Fcor.put("o1",oct); Fcor.put('d1',dct)
Fcor.put("n2",ncx); Fcor.put("o2",ocx); Fcor.put('d2',dcx)
Fcor.put("n3",ncy); Fcor.put("o3",ocy); Fcor.put('d3',dcy)    
Fcor.put("n4",ncz); Fcor.put("o4",ocz); Fcor.put('d4',dcz)    

cor = np.zeros(2*nt-1,'f')
       
# ------------------------------------------------------------

for icz in range(-nzlag,+nzlag+1): 
    for icy in range(-nylag,+nylag+1):
        for icx in range(-nxlag,+nxlag+1):
        
            cor = np.correlate(wfl[izmid-icz,iymid-icy,izmid-icz,:],
                               opr[izmid+icz,iymid+icy,ixmid+icx,:],mode='full')
            Fcor.write(cor[nt-ntlag:nt+ntlag+1])

# ------------------------------------------------------------
Fopr.close()
Fwfl.close()
Fcor.close()

