#!/usr/bin/env python
import sys, math
import vplot

def arrow(vp,x1,y1,x2,y2):
    dx = x2 - x1
    dy = y2 - y1
    r = math.hypot(dx,dy)
    if  r < .5:  
        r = .5

    backx = -.15 * dx / r
    backy = -.15 * dy / r
    perpx =  .05 * dy / r
    perpy = -.05 * dx / r
    vp.umove( x1, y1 )         
    vp.udraw( x2, y2 )

    tipx = x2 + backx + perpx
    tipy = y2 + backy + perpy
    vp.umove( x2, y2 )         
    vp.udraw( tipx, tipy)

    tipx = x2 + backx - perpx
    tipy = y2 + backy - perpy
    vp.umove( x2, y2 )         
    vp.udraw( tipx, tipy)

def filter(vp,x0,y0,letter):
    dx=0.5
    dy=0.5
    rad=2.0
    e=0.4

    # vertical half line
    vp.move( x0+.5*dx,  y0-.5*dy)
    vp.draw( x0+.5*dx,  y0+(2.5+e)*dy)

    # vertical lines
    vp.move( x0- .5*dx, y0-(2.5+e)*dy)
    vp.draw( x0- .5*dx, y0+(2.5+e)*dy)

    vp.move( x0-1.5*dx, y0-(2.5+e)*dy)
    vp.draw( x0-1.5*dx, y0+(2.5+e)*dy)

    vp.move( x0-2.5*dx, y0-(2.5+e)*dy)
    vp.draw( x0-2.5*dx, y0+(2.5+e)*dy)

    # short horiz lines
    vp.move( x0-(2.5+e)*dx, y0-2.5*dy)
    vp.draw( x0-     .5*dx, y0-2.5*dy)

    vp.move( x0-(2.5+e)*dx, y0-1.5*dy)
    vp.draw( x0-     .5*dx, y0-1.5*dy)

    # longer horiz lines
    vp.move( x0-(2.5+e)*dx, y0- .5*dy)
    vp.draw( x0+     .5*dx, y0- .5*dy)

    vp.move( x0-(2.5+e)*dx, y0+ .5*dy)
    vp.draw( x0+     .5*dx, y0+ .5*dy)

    vp.move( x0-(2.5+e)*dx, y0+1.5*dy)
    vp.draw( x0+     .5*dx, y0+1.5*dy)

    vp.move( x0-(2.5+e)*dx, y0+2.5*dy)
    vp.draw( x0+     .5*dx, y0+2.5*dy)

    vp.penup()

    # circles
    for theta in range(90,280,10):
        x = x0 + rad * math.cos( theta * math.pi/180)
	y = y0 + rad * math.sin( theta * math.pi/180)
        vp.pendn(x,y)

    # arrows
    for theta in range(90+180/10,280,180/5):
        x  = x0 + rad      * math.cos( theta * math.pi/180)
        y  = y0 + rad      * math.sin( theta * math.pi/180)
        xp = x0 + (rad+.5) * math.cos( theta * math.pi/180)
        yp = y0 + (rad+.5) * math.sin( theta * math.pi/180)
        arrow(vp,x,y,xp,yp)

    vp.text(x0-.3*dx, y0-.3*dy, 10, 0, letter)

if __name__ == "__main__":
    vp = vplot.Vplot()
    vp.fat(3)
    
    vp.color(7)
    filter(vp,8,4,'A')
    
    vp.color(6)
    filter(vp,6,6,'B')

    sys.exit(0)



