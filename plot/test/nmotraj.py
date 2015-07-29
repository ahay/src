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
import math
import vplot

def onebox(vp, x, y, dx, dy):
    edge=0.95
    vx = []
    vy = []

    vx.append(x-edge*dx) 
    vy.append(y-edge*dy)
    
    vx.append(x+edge*dx)      
    vy.append(y-edge*dy)
    
    vx.append(x+edge*dx)      
    vy.append(y+edge*dy)
    
    vx.append(x-edge*dx)      
    vy.append(y+edge*dy)

    vp.umove( vx[3], vy[3])
    vp.udraw( vx[0], vy[0])
    vp.udraw( vx[1], vy[1])
    vp.udraw( vx[2], vy[2])
    vp.udraw( vx[3], vy[3])

def draw(vp):
    nt=48
    nx=18
    xmin=-2
    xmax=9
    tmax=9.5
    
    dt =  tmax / (nt-1.0)
    dx =  xmax / (nx-1.0)
    t0 = 0
    x0 = dx/2
    z = 0
    v = 1.01 * xmax/math.sqrt( tmax*tmax - z*z)

    vp.uorig( -1.+xmin, 0.)
    vp.umove( xmin, 9.5)       
    vp.udraw( xmax, 9.5)
    vp.umove(   0., 9.5)     
    vp.udraw( 0., 9.5-tmax)

    vp.color(5)
    vp.utext( xmax-.35, 9.5-.45, 12, 0, "x")
    vp.utext( .25     ,  .2    , 12, 0, "t")

    vp.color(6)
    for iz in range(1,4):
        z = (-.10 + .33*iz) * tmax + dt/2
        itlast = int(z / dt)
        
        for ix in xrange(nx):
	    x = x0 + ix*dx
	    t = math.hypot(z,x/v)
	    it = int(t / dt)
	    t = t0+it*dt

            if it< itlast+1:
                i0=it
            else:
                i0=itlast+1

        
            for i in xrange(i0,it+1):
      		t = t0+i*dt
                if t < tmax:
		    onebox(vp,  x, 9.5-t, dx/2, dt/2 )
                    if -x > xmin:
                        onebox(vp, -x, 9.5-t, dx/2, dt/2 )
            itlast = it

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp)

