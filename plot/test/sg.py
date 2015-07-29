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
import sys
import vplot

def draw(vp,dgods):
    s0=0.5
    g0=0.5

    ns=13
    size=10
    nh=10
	
    ds= (10.24-1.0)/ns
    dg = ds * dgods
    vp.clip( 0.,0., 10.24/.75, 10.24)

    vp.move( g0, 1.)
    vp.draw( g0, 8.5)
    vp.text( g0, 9.,  size+4, 0, "s")

    vp.move( 8., s0)
    vp.draw( 12., s0)
    vp.text( 12.5, s0, size+4, 0, "g")

    for i in range(ns):
	s = s0 + ds*i
	for ih in range(nh):  
	    h = dg*ih
	    g = s + h
	    if h < 0.1+ 9*ds:
		vp.text( g, s, size, 0, str(ih))

if __name__ == "__main__":
    vp = vplot.Vplot()
    dgods = float(sys.argv[1])
    draw(vp,dgods)
