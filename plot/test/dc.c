/* */
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

static void doit(float wide, int panel);

int main (void)
{
    int j;
    float wide=1.4;

    vp_init();

    for (j=0; j <= 3; j++) {
	vp_uorig ( -.5-wide - 2.3*wide * j, 0.);
	doit( wide,j);
    }

    return 0;
}

static void doit(float wide, int panel)
{
    float v, t,z,x,H, tmax=4., dx=.04;
    int	labelsize=7, iz;

    H=7.; /* height of top of frame */
    v= 1.3 * wide/tmax;

    vp_fat (6);
    vp_umove (-wide,H-0.); 	
    vp_udraw (wide,H-0);
    vp_umove (0,H-0.);		
    vp_udraw (0.,H-tmax);

    vp_fat (3);
    vp_utext( wide-.3, H-       0.4, labelsize, 0, "x");	
    vp_utext(     0.15,H- tmax+.03 , labelsize, 0, "t");

    switch (panel) {
	case 0:
	    vp_utext( -wide+.1, H+       0.4, labelsize, 0, "at z0 (beach)");
	    break;
	case 1:
	    vp_utext( -wide+.1, H+       0.4, labelsize, 0, "at z1");
	    break;
	case 2:
	    vp_utext( -wide+.1, H+       0.4, labelsize, 0, "at z2");
	    break;
	case 3:
	    vp_utext( -wide+.1, H+       0.4, labelsize, 0, "at z3 (barrier)");
	    break;
    }

    iz = 3 - panel;
    z    = .22 * iz * tmax;
    vp_penup ();
    for( x=-wide+dx/2.; x<wide; x += dx) {
	t = hypotf(z,x/v);
	if( t< sqrtf(2.)*z && t< tmax)
	    vp_upendn (x, H-t);
    }
}

/* 	$Id: dc.c 7107 2011-04-10 02:04:14Z ivlad $	 */
