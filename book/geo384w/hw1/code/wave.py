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
        self.c0  = -2.0*(self.c11+self.c12+self.c21+self.c22)
    
    def apply(self,uin,uout):
        n1,n2 = uin.shape
    
        uout[2:n1-2,2:n2-2] = \
         self.c11*(uin[1:n1-3,2:n2-2]+uin[3:n1-1,2:n2-2]) + \
         self.c12*(uin[0:n1-4,2:n2-2]+uin[4:n1  ,2:n2-2]) + \
         self.c21*(uin[2:n1-2,1:n2-3]+uin[2:n1-2,3:n2-1]) + \
         self.c22*(uin[2:n1-2,0:n2-4]+uin[2:n1-2,4:n2  ]) + \
         self.c0*uin[2:n1-2,2:n2-2]

par = m8r.Par()

# setup I/O files
Fr=m8r.Input()       # source position
Fo=m8r.Output()      # output wavefield

Fv=m8r.Input ("v")   # velocity
Fw=m8r.Input ("wav") # source wavefield


# Read/Write axes
a1 = Fr.axis(1); n1 = a1['n']; d1 = a1['d']
a2 = Fr.axis(2); n2 = a2['n']; d2 = a2['d']
at = Fw.axis(1); nt = at['n']; dt = at['d']

ft = par.int('ft',0)
jt = par.int('jt',0)

Fo.put('n3',(nt-ft)/jt)
Fo.put('d3',jt*dt)
Fo.put('o3',ft*dt)

dt2 =    dt*dt

# set Laplacian coefficients 
laplace = Laplacian(1.0/(d1*d1),1.0/(d2*d2))

# read wavelet, velocity & source position
ww = numpy.zeros(nt,'f');      Fw.read(ww)
vv = numpy.zeros([n2,n1],'f'); Fv.read(vv)
rr = numpy.zeros([n2,n1],'f'); Fr.read(rr)

# allocate temporary arrays
u0 = numpy.zeros([n2,n1],'f')
u1 = numpy.zeros([n2,n1],'f')
u2 = numpy.zeros([n2,n1],'f')
ud = numpy.zeros([n2,n1],'f')

vv = vv*vv*dt2

# Time loop
for it in range(nt):
    laplace.apply(u1,ud)

    # scale by velocity 
    ud = ud*vv
    # inject wavelet 
    ud = ud + ww[it] * rr
    # time step 
    u2 = 2*u1 - u0 + ud 
    u0 = u1
    u1 = u2

    # write wavefield to output 
    if it >= ft and 0 == (it-ft)%jt:
        sys.stderr.write("\b\b\b\b\b %d" % it)
        Fo.write(u1)
