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

def arrow(vp,x1,y1,x2,y2):
    n=200

    dx = x2 - x1
    dy = y2 - y1
    r = math.hypot(dx,dy)
    if r < 0.5:
        r = 0.5

    backx = -0.40 * dx / r
    backy = -0.40 * dy / r
    perpx =  0.14 * dy / r
    perpy = -0.14 * dx / r

    vp.umove(x1,y1)
    noise = 1/600.
    for i in xrange(n):
	vp.udraw(x1 + ((i+1)*(x2-x1))/n + noise, 
		 y1 + ((i+1)*(y2-y1))/n + noise)
	noise = -noise

    tipx = x2 + backx + perpx
    tipy = y2 + backy + perpy
    vp.umove(x2, y2)		
    vp.udraw(tipx,tipy)

    tipx = x2 + backx - perpx
    tipy = y2 + backy - perpy
    vp.umove(x2, y2)		
    vp.udraw(tipx, tipy)

def draw(vp):
    top=10.
    xmin=0.
    ymin=0.
    xmax = 13.
    ymax = 3.
    d = 2.

    vp.uorig (-.1,-.1)
    vp.color (4)

    # axes 
    vp.fat(7)
    vp.umove( xmin, top-  0.)  
    vp.udraw( xmax, top-  0.)
    vp.umove( 0.,   top-ymin)  
    vp.udraw( 0.,   top-ymax)
    vp.fat(3)

    # rays 
    vp.color(7)
    vp.umove( xmin,           top-ymin)
    vp.udraw( xmin+d,         top-( ymin+d ))	# diagonal down 
    vp.udraw( xmin+d+d+d+d+d, top-( ymin+d ))	# along bed 

    arrow(vp,
          xmin+3*d, top-(ymin+d),
	  xmin+4*d, top-(ymin  ))
    arrow(vp,
          xmin+4*d, top-(ymin+d),
	  xmin+5*d, top-(ymin  ))
    arrow(vp,
          xmin+5*d, top-(ymin+d),
	  xmin+6*d, top-(ymin  ))

if __name__ == "__main__":
    vp = vplot.Vplot()
    draw(vp)

