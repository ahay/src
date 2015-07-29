#!/usr/bin/env python
##   Copyright (C) 2008 University of Texas at Austin
##  
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##  
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##  
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
import math, sys
import numpy
import m8r, vplot

def main(argv=sys.argv):

    nt=300
    nx=64
    nz=300
    nb=50
    v=1.0
    top=3.2
    c1=0.9
    c2=6.8
    d1=0.02
    d2=0.12
    
    par = m8r.Par(argv)    

    top = par.float('top',5.0)
    c1 = par.float('c1',0.5)
    c2 = par.float('c2',5.0)

    vp = vplot.Vplot()

    vp.uorig (-c1,-.5)
    vp.uclip (0.,top-3.,4.,top)
    vp.umove (0.,top-0.)  
    vp.udraw (0.,top-4.)
    vp.umove (0.,top-0.)  
    vp.udraw (4.,top-0.)

    z = 0.4
    while z < 4: 
	vp.penup()
	x0 = z * math.tan(math.pi*45./180.)
        x = 0
        while x < 4:
	    t = math.hypot(z,x-x0)/v
	    vp.upendn (x,top-t)
            x += 0.01
        z += 0.4

    b = numpy.zeros(nb,'f')
    for ib in xrange(nb):
	b[ib] = math.exp(-3.*(ib+1.)/20.) * math.sin(math.pi*(ib+1)/10)

    cs = par.string('c')
    if cs:
	c = m8r.Output('c')
	c.setformat('native_float')
	c.put('n1',nt)
	c.put('n2',nx)
	c.put('d1',d1)
	c.put('d2',d2)
	
        tdat = numpy.zeros((nx,nt),'f')

        for iz in range(12):
	    z = (iz+1.)*nt/12.
	    x0 = z * math.tan( math.pi*45./180.)
            for ix in xrange(nx):
		x = (ix+1.) * d2/d1
		t = math.hypot(z,x-x0)/v
                for ib in xrange(nb):
		    it = t+ib-1
		    if it < nt:
                        tdat[ix,it] += b[ib]
        c.write(tdat)

    vp.uorig (-c2,-.5)
    vp.uclip (0.,top-3.,4.,top)
    vp.umove (0.,top-0.)  
    vp.udraw (0.,top-4.)
    vp.umove (0.,top-0.)   
    vp.udraw (4.,top-0.)

    t = 0.4
    while t < 6:
	vp.penup ()
	x0 = t / math.sin( math.pi*45./180.)
        theta=-89.5
	while theta<89.5:
	    z =      t * math.cos (math.pi*theta/180)
	    x = x0 + t * math.sin (math.pi*theta/180)
	    vp.upendn (x,top-z)
            theta += 1
        t += 0.4

    ds = par.string('d')
    if ds:
	d = m8r.Output('d')
	d.setformat('native_float')
	d.put('n1',nz)
	d.put('n2',nx)
	d.put('d1',d1)
	d.put('d1',d2)
	
        zdat = numpy.zeros([nx,nz],'f')

	for it in range(20):
	    t = (it+1.)*nz/20.
	    x0 = t / math.sin( math.pi*45./180.)
            theta=-89.5
            while theta<89.5:
		z =      t * math.cos (math.pi*theta/180.)
		x = x0 + t * math.sin (math.pi*theta/180.)
		ix = x * d1/d2
		r = math.hypot(z,x-x0)
                for ib in xrange(nb):
		    iz = z + ib-1
		    if iz >= 0 and iz < nz and ix >=0 and ix < nx:
			zdat[ix,iz] += b[ib]*r
                theta += 1
	
        d.write(zdat)

if __name__ == '__main__':
    main()
