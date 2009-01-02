#!/usr/bin/env python
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
