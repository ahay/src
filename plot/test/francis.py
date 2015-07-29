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
import math, cmath
import vplot

def draw(vp,nfat=1,ntheta=271,eps=0.1,k2=0.05):
    for job in range(5):
	vp.uorig ( -1.75-2.20*job, -1.5 )
	vp.uclip (-1.,-1.,1.3,1.1)
	vp.fat (nfat)

	vp.umove (0.,-.9)  
	vp.udraw (0.,4.9)
	vp.umove (-.5,0.)  
	vp.udraw (1.2,0.)

	vp.utext(  .1,  .9,  4,0, "Im")
	vp.utext( 1.05, -.2, 4,0, "Re")

	if job == 0:
            vp.utext( .1, -.7 , 4,0, "\\F10 w=+p")
            vp.utext( .1,  .6 , 4,0, "\\F10 w=-p")
            vp.utext( .1,  .05, 4,0, "\\F10 w=0")
        elif job == 1:
            vp.utext( .4, -.65, 4,0, "\\F10 w=p/2")
            vp.utext( .4,  .55, 4,0, "\\F10 w=-p/2")
            vp.utext( .1,  .05, 4,0, "\\F10 w=0")


	vp.fat(0)
	vp.penup()

        theta = -180
        while theta<180.1:
	    rads = 2. * math.pi * theta / 360.
	    cz = cmath.exp(rads*1j)
	    cs = (1+eps - cz )/2

            cs = (lambda c: 0.05+0.8*theta*1j/180,
                  lambda c: c,
                  lambda c: c*c,
                  lambda c: c*c + k2,
                  lambda c: cmath.sqrt(c*c + k2))[job](cs)

	    x = cs.real 
	    y = cs.imag

	    vp.upendn ( x, -y)
            theta += 360.0/ntheta

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp)


