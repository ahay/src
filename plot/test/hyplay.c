/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>

#include <rsfplot.h>

int main (void)
{
    float t,x,y;
    const float v=1., dx=.02, dt=.2, tmax=4., dy=.3, H=6., z=2.5;
    
    vp_init();
    vp_uorig ( -3., 0.);

    vp_fat ( 2);
    vp_umove (-2.,H);  
    vp_udraw (+2.,H);
    vp_umove (0,H);  
    vp_udraw (0.,H-tmax);

    vp_fat ( 0);
    vp_utext( 1.8, H-       0.3, 5, 0, "x");
    vp_utext( 0.15,H- tmax+.03 , 5, 0, "t");

    vp_uclip ( -2.-dx/2 ,H-tmax-dt/2+.5, 2.+dx/2, H+.05);

    for(y=-4.+dy/2.; y<4.; ) { 
	vp_penup();
	for (x=-2.+dx/2.; x < 2.; x += dx, y += dy) { 
	    t = hypotf(z,x-y)/v;
	    vp_upendn (x, H-t);
	}
    }

    return 0;
}

