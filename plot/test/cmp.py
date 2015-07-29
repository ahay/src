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

iseed=1996

def rand01(ia = 727,im = 524287):
    global iseed
    iseed = (iseed*ia)%im
    return (float(iseed) - 0.5)/(float(im) - 1.0)

def draw(vp):
    xmin=-4.0
    xmax=7.0
    tmax=9.0
    v=1.0
    nz=100
    nx=100
    degrees=85

    t0 = 0	
    x0= xmin
    dx= (xmax-xmin)/(nx-1)

    vp.uorig( -1.+xmin, 0.)

    vp.utext(  xmin/2, tmax+.45, 18, 0, "Common Midpoint")
    vp.utext( xmax+.1, tmax-.45, 14, 0, "h")
    vp.utext( .25     ,  .1    , 14, 0, "t")

    vp.uclip (xmin, .6, xmax, tmax)
    vp.umove( xmin, tmax) 	
    vp.udraw( xmax, tmax)
    vp.umove(   0., tmax) 	
    vp.udraw( 0., tmax-tmax)

    vp.color(6)
    for iz in xrange(nz):
	alfa = degrees *  math.pi / 180 * rand01() 
	ca = math.cos( alfa)
	sa = math.sin( alfa)
	z =  tmax * (.1 + .9*rand01())
	z =  tmax * rand01()
	vp.penup()
        for ix in xrange(nx):
	    x = x0 + ix*dx
	    t = math.hypot( z,ca*x/v)
	    if x > -v*t and x <  v*t:
		vp.upendn(x, tmax-t)

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp)
