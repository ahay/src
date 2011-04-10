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

static void doit(float wide, int job);

int main (void)
{
    float wide=3.;

    vp_init();

    vp_uorig ( -.5-wide, 0.);		
    doit( wide, 1 );
    vp_uorig ( -.5-wide-2.2*wide, 0.);	
    doit( wide, 2 );

    return 0;
}

static void doit(float wide, int job)
{
    float v, t,z, x,H, tmax=8., zbot,dx=.04;
    int	labelsize=12, iz;

    H=9.;			/* height of top of frame */

    v= 1.3 * wide/tmax;

    vp_fat ( 6);
    vp_umove (-wide,H-0.); 	
    vp_udraw (wide,H-0);
    vp_umove (0,H-0.);		
    vp_udraw (0.,H-tmax);

    vp_fat ( 3);
    vp_utext( wide-.3, H-       0.4, labelsize, 0, "x");	
    if (job == 1) {
	vp_utext(     0.15,H- tmax+.03 , labelsize, 0, "t");	
    } else {
	vp_utext(     0.15,H- tmax+.03 , labelsize, 0, "t+z/v");
    }

    for (iz=1; iz <= 3; iz++) {
	z   = .2*iz*tmax;
	zbot= .2* 3*tmax;
	vp_penup ();
	for( x=-wide+dx/2.; x<wide; x += dx) {
	    t = hypotf(z,x/v);
	    if( t < sqrtf(2.)*z) {
		if( job == 2) t += zbot-z;
		vp_upendn (x, H-t);
	    }
	}
    }
}

/* 	$Id$	 */
