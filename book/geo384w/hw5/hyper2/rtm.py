#!/usr/bin/env python

import sys
import numpy
import m8r

class Laplacian:
    'Laplacian operator, 4th-order finite-difference'
    
    def __init__(self,d1,d2):
        self.c11 = 4.0*d1/3.0
        self.c12=  -d1/12.0
        self.c21 = 4.0*d2/3.0
        self.c22=  -d2/12.0
        self.c1  = -2.0*(self.c11+self.c12)
        self.c2  = -2.0*(self.c21+self.c22)
    
    def apply(self,uin,uout):
        n1,n2 = uin.shape
    
        uout[2:n1-2,2:n2-2] = \
         self.c11*(uin[1:n1-3,2:n2-2]+uin[3:n1-1,2:n2-2]) + \
         self.c12*(uin[0:n1-4,2:n2-2]+uin[4:n1  ,2:n2-2]) + \
         self.c1*uin[2:n1-2,2:n2-2] + \
         self.c21*(uin[2:n1-2,1:n2-3]+uin[2:n1-2,3:n2-1]) + \
         self.c22*(uin[2:n1-2,0:n2-4]+uin[2:n1-2,4:n2  ]) + \
         self.c2*uin[2:n1-2,2:n2-2]

par = m8r.Par()

# setup I/O files
modl=m8r.Input()       # velocity model
imag=m8r.Output()      # output image

data=m8r.Input ('data') # seismic data
wave=m8r.Output('wave') # wavefield

# Dimensions
n1 = modl.int('n1')
n2 = modl.int('n2')

dz = modl.float('d1')
dx = modl.float('d2')

nt = data.int('n1')
dt = data.float('d1')

nx = data.int('n2')
if nx != n2:
    raise RuntimeError('Need n2=%d in data',n2)

n0 = par.int('n0',0) # surface
jt = par.int('jt',1) # time interval

wave.put('n3',1+(nt-1)/jt)
wave.put('d3',-jt*dt)
wave.put('o3',(nt-1)*dt)

dt2 =    dt*dt

# set Laplacian coefficients 
laplace = Laplacian(1.0/(dz*dz),1.0/(dx*dx))

# read data and velocity
dd = numpy.zeros([n2,nt],'f')
data.read(dd)
vv = numpy.zeros([n2,n1],'f')
modl.read(vv)

# allocate temporary arrays
u0 = numpy.zeros([n2,n1],'f')
u1 = numpy.zeros([n2,n1],'f')
ud = numpy.zeros([n2,n1],'f')

vv *= vv*dt2

# Time loop
for it in range(nt-1,-1,-1):
    sys.stderr.write("\b\b\b\b\b %d" % it)
    
    laplace.apply(u1,ud)

    # scale by velocity
    ud *= vv

    # time step 
    u2 = 2*u1 - u0 + ud 
    u0 = u1
    u1 = u2

    # inject data
    u1[:,n0] += dd[:,it]

    if 0 == it%jt:
        wave.write(u1)

# output image
imag.write(u1)
