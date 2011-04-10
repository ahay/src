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

static void onebox(float x, float y, float dx, float dy );

int main (void)
{
    int it,nt=48, ix,nx=18, i, itlast, iz;
    float dt,dx, xmin=-2., xmax=9., tmax=9.5, x,t,z, x0, t0, v;

    dt =  tmax / (nt-1);
    dx =  xmax / (nx-1);
    t0 = 0.;
    x0 = dx/2;
    z = 0.; /* ? */
    v = 1.01 * xmax/ sqrtf( tmax*tmax - z*z);

    vp_init();

    vp_uorig( -1.+xmin, 0.);
    vp_umove( xmin, 9.5);       
    vp_udraw( xmax, 9.5);
    vp_umove(   0., 9.5);       
    vp_udraw( 0., 9.5-tmax);

    vp_color(5);
    vp_utext( xmax-.35, 9.5-.45, 12, 0, "x");
    vp_utext( .25     ,  .2    , 12, 0, "t");

    vp_color(6);
    for (iz=1; iz <= 3; iz++) {
        z = (-.10 + .33*iz) * tmax + dt/2;
        itlast = z / dt;
        for (ix=0; ix < nx; ix++) {
	    x = x0 + ix*dx;
	    t = hypotf(z,x/v);
	    it = t / dt;
	    t = t0+it*dt;

	    for (i = (it< itlast+1)? it: itlast+1; i <= it; i++) {
		t = t0+i*dt;
		if( t < tmax ) {
		    onebox(  x, 9.5-t, dx/2, dt/2 );
		    if( -x > xmin )
			onebox( -x, 9.5-t, dx/2, dt/2 );
		}
	    }
	    itlast = it;
	}
    }

    return 0;
}

static void onebox(float x, float y, float dx, float dy )
{
    float vx[4],vy[4];
    const float edge=0.95;

    vx[0] = x-edge*dx; 
    vy[0] = y-edge*dy;
    
    vx[1] = x+edge*dx;      
    vy[1] = y-edge*dy;
    
    vx[2] = x+edge*dx;      
    vy[2] = y+edge*dy;
    
    vx[3] = x-edge*dx;      
    vy[3] = y+edge*dy;

    vp_umove( vx[3], vy[3]);
    vp_udraw( vx[0], vy[0]);
    vp_udraw( vx[1], vy[1]);
    vp_udraw( vx[2], vy[2]);
    vp_udraw( vx[3], vy[3]);
}

/* 	$Id$	 */
