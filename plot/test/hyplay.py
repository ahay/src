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

def draw(vp):
    v=1
    dx=0.02
    dt=0.2
    tmax=4
    dy=0.3
    H=6
    z=2.5
    
    vp.uorig( -3., 0.)

    vp.fat( 2)
    vp.umove(-2.,H)  
    vp.udraw(+2.,H)
    vp.umove(0,H)  
    vp.udraw(0.,H-tmax)

    vp.fat(0)
    vp.utext( 1.8, H-       0.3, 5, 0, "x")
    vp.utext( 0.15,H- tmax+.03 , 5, 0, "t")

    vp.uclip ( -2.-dx/2 ,H-tmax-dt/2+.5, 2.+dx/2, H+.05)

    y=-4.+dy/2.
    while y < 4:
	vp.penup()
        x=-2.+dx/2.
        while x < 2:
            t = math.hypot(z,x-y)/v
	    vp.upendn (x, H-t)
            x += dx
            y += dy

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp)




