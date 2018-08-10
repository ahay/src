#!/usr/bin/env python

import sys
import numpy
import m8r

c0=-30./12.
c1=+16./12.
c2=- 1./12.

par = m8r.Par()
verb = par.bool("verb",False) # verbosity

# setup I/O files
Fw=m8r.Input()
Fv=m8r.Input ("vel")
Fr=m8r.Input ("ref")
Fo=m8r.Output()

# Read/Write axes
at = Fw.axis(1); nt = at['n']; dt = at['d']
az = Fv.axis(1); nz = az['n']; dz = az['d']
ax = Fv.axis(2); nx = ax['n']; dx = ax['d']

Fo.putaxis(az,1)
Fo.putaxis(ax,2)
Fo.putaxis(at,3)

dt2 =    dt*dt
idz = 1/(dz*dz)
idx = 1/(dx*dx) 

# read wavelet, velocity & reflectivity
ww = numpy.zeros(nt,'f');      Fw.read(ww)
vv = numpy.zeros([nz,nx],'f'); Fv.read(vv)
rr = numpy.zeros([nz,nx],'f'); Fr.read(rr)

# allocate temporary arrays
um = numpy.zeros([nz,nx],'f')
uo = numpy.zeros([nz,nx],'f')
up = numpy.zeros([nz,nx],'f')
ud = numpy.zeros([nz,nx],'f')

# MAIN LOOP
for it in range(nt):
    if verb:
        sys.stderr.write("\b\b\b\b\b %d" % it)

    ud[2:-2,2:-2] = \
    c0* uo[2:-2,2:-2] * (idx + idz)        + \
    c1*(uo[2:-2,1:-3] + uo[2:-2,3:-1])*idx + \
    c2*(uo[2:-2, :-4] + uo[2:-2,4:  ])*idx + \
    c1*(uo[1:-3,2:-2] + uo[3:-1,2:-2])*idz + \
    c2*(uo[ :-4,2:-2] + uo[4:  ,2:-2])*idz

    # inject wavelet
    ud = ud - ww[it] * rr

    # scale by velocity
    ud= ud *vv*vv

    # time step
    up = 2*uo - um + ud * dt2
    um =   uo
    uo =   up

    Fo.write(uo)
if verb:
    sys.stderr.write("\n")

sys.exit(0)
