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
    nx=100
    nz=100
    degrees=85
    xmin=-4.0
    xmax=7.0
    tmax=9

    x0= xmin
    dx= (xmax-xmin)/(nx-1)
    vp.uorig( -1.+xmin, 0.)

    vp.utext(  xmin/2, tmax+.45, 18, 0, "Common Shot")
    vp.utext( xmax+.1, tmax-.45, 14, 0, "g")
    vp.utext( .25     ,  .1    , 14, 0, "t")

    vp.uclip (xmin, .6, xmax, tmax)
    vp.umove( xmin, tmax) 	
    vp.udraw( xmax, tmax)
    vp.umove(   0., tmax) 	
    vp.udraw( 0., tmax-tmax)
    vp.color(6)
    for iz in xrange(nz):
	orig = 3 * ((xmax-xmin) * rand01() +xmin)
	alfa = degrees *  2 * 3.14 / 360 * rand01() 
	c2a = math.cos( 2*alfa)
	vp.penup()
        for ix in xrange(nx):       # x=g s=0 
	    x = x0 + ix*dx
	    arg = orig*orig +(x-orig)*(x-orig) + 2*orig*(x-orig) * c2a
	    t = math.sqrt( arg) 
	    if t < math.fabs(x):
                t = math.fabs(x)
            vp.upendn(x, tmax-t)

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp)
