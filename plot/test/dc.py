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

def doit(vp,wide,panel):
    tmax=4.
    dx=.04
    labelsize=7

    H=7. # height of top of frame 
    v= 1.3 * wide/tmax

    vp.fat (6)
    vp.umove (-wide,H-0.) 	
    vp.udraw (wide,H-0)
    vp.umove (0,H-0.)		
    vp.udraw (0.,H-tmax)

    vp.fat (3)
    vp.utext( wide-.3, H-       0.4, labelsize, 0, "x")	
    vp.utext(     0.15,H- tmax+.03 , labelsize, 0, "t")

    (lambda: vp.utext( -wide+.1, H+       0.4, labelsize, 0, "at z0 (beach)"),
     lambda: vp.utext( -wide+.1, H+       0.4, labelsize, 0, "at z1"),
     lambda: vp.utext( -wide+.1, H+       0.4, labelsize, 0, "at z2"),
     lambda: vp.utext( -wide+.1, H+       0.4, labelsize, 0, "at z3 (barrier)")
     )[panel]()

    iz = 3 - panel
    z    = .22 * iz * tmax
    vp.penup ()
    x=-wide+dx/2
    while (x<wide):
	t = math.hypot(z,x/v)
	if t< math.sqrt(2)*z and t< tmax:
	    vp.upendn (x, H-t)
        x += dx

def draw(vp,wide=1.4):
    for j in range(4):
	vp.uorig ( -.5-wide - 2.3*wide * j, 0.)
	doit(vp, wide,j)

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp)
