/* Simple v(z) synthetic.

Notes about theory:

v = v0 + alpha * z
t = 2 * \int_0^z dz/(v0+alpha*z)
t = 2 * log( 1 + alpha*z/v0 ) / alpha
exp(alpha*tmax/2.) = 1 + alpha * zmax/ v0
v0 * (exp( alpha * tmax/2.) - 1) /alpha =  zmax = dz * (nz+1)
*/
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

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt, nx, it,ix;
    float t,t0,dt, x,x0,dx, z; 
    float tmax,xmax,delx, dxdz,v0,alpha;
    float** datr;
    sf_file sag;

    sf_init (argc, argv);
    sag = sf_output("out");

    if (!sf_getint ("nt",&nt)) nt = 200;
    /* Number of samples in time */
    if (!sf_getint ("nx",&nx)) nx = 200;
    /* Number of samples in distance */

    if (!sf_getfloat ("tmax",&tmax)) tmax = 4.;
    /* Maximum time */
    if (!sf_getfloat ("xmax",&xmax)) xmax = 4.;
    /* Maximum distance */
    if (!sf_getfloat ("delx",&delx)) delx = .5;
    /* Increment in x */
    if (!sf_getfloat ("dxdz",&dxdz)) dxdz = 1.;
    /* Slope for the line of diffractors */
    if (!sf_getfloat ("v0",&v0)) v0 = 1.5;
    /* Initial velocity */
    if (!sf_getfloat ("alpha",&alpha)) alpha = 0.;
    /* Velocity gradient */

    t0 = 0;		
    x0 = 0.;
    dt = tmax / nt;
    dx = xmax / nx;

    sf_setformat(sag,"native_float");    
    sf_putint (sag,"n1",nt); 
    sf_putfloat (sag,"o1",t0); 
    sf_putfloat (sag,"d1",dt);
    sf_putint (sag,"n2",nx); 
    sf_putfloat (sag,"o2",x0); 
    sf_putfloat (sag,"d2",dx);
    sf_putstring (sag,"label1","Pseudo-depth (s)");
    sf_putstring (sag,"label2","Lateral (km)");

    datr = sf_floatalloc2(nt,nx);
    for (ix=0; ix < nx; ix++) {
	for (it=0; it < nt; it++) {
	    datr[ix][it] = 0.;
	}
    }

    for (x = delx/2; x <= xmax; x+= delx) {
	z = x / dxdz;
	if( alpha != 0.) {
	    t = 2. * log( 1. + alpha * z / v0) / alpha;
	} else {		 
	    t = 2. * z / v0;
	}
	it = 0.5 + t/dt;
	ix = 0.5 + x/dx;
	if(ix < nx && it < nt) datr[ix][it] += 1.;
    }

    sf_floatwrite (datr[0],nt*nx,sag);

    exit (0);
}
