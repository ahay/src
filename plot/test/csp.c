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

static int iseed=1996;

static float rand01(void);

int main(void)
{
    int ix,nx=100, iz, nz=100, degrees=85;
    float dx, xmin=-4., xmax=7., tmax=9., x,t, x0;
    float alfa, c2a, orig, arg;

    vp_init();

    x0= xmin;
    dx= (xmax-xmin)/(nx-1);
    vp_uorig( -1.+xmin, 0.);

    vp_utext(  xmin/2, tmax+.45, 18, 0, "Common Shot");
    vp_utext( xmax+.1, tmax-.45, 14, 0, "g");
    vp_utext( .25     ,  .1    , 14, 0, "t");

    vp_uclip (xmin, .6, xmax, tmax);
    vp_umove( xmin, tmax); 	
    vp_udraw( xmax, tmax);
    vp_umove(   0., tmax); 	
    vp_udraw( 0., tmax-tmax);
    vp_color(6);
    for (iz=0; iz < nz; iz++) {
	orig = 3 * ((xmax-xmin) * rand01() +xmin);
	alfa = degrees *  2 * 3.14 / 360 * rand01(); 
	c2a = cos( 2*alfa);
	vp_penup();
	for (ix=0; ix < nx; ix++) {			/* x=g; s=0 */
	    x = x0 + ix*dx;
	    arg = orig*orig +(x-orig)*(x-orig) + 2*orig*(x-orig) * c2a;
	    t = sqrtf( arg); 
	    if( t < fabsf(x)) t = fabsf(x);
	    vp_upendn(x, tmax-t);
	}
    }

    return 0;
}

float rand01(void)
{
    const int ia = 727, im = 524287;
    iseed = (iseed*ia)%im;
    return ((float) iseed - 0.5)/((float) im - 1.);
}

/* 	$Id: csp.c 7107 2011-04-10 02:04:14Z ivlad $	 */
