/* Conversion from depth to time in a V(z) medium.

Transforms function of z to function of

tau = Integral[2/v[x,n],{n,0,z}]
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

#include <rsf.h>

int main (int argc, char *argv[])
{
    int nt, nz, nx, iz, ix;
    bool slow;
    float t0, dt, z0, dz, z=0., eps;
    float *time, *depth, *vel;
    sf_map str;
    sf_file in, out, velocity;

    sf_init(argc, argv);
    in       = sf_input("in");
    velocity = sf_input("velocity");
    out      = sf_output("out");

    if (!sf_histint (in,"n1",&nz)) sf_error ("No n1= in input");
    if (!sf_histfloat (in,"d1",&dz)) sf_error ("No d1= in input");
    if (!sf_histfloat (in,"o1",&z0)) z0 = 0.;

    nx = sf_leftsize(in,1);
    
    if (!sf_getint ("nt",&nt)) {
	/* Number of points in time (default is n1) */
	nt = nz; 
    } else {
	sf_putfloat(out,"n1",nt);
    }
    if (!sf_getfloat ("dt",&dt)) {
	/* Time sampling (default is d1) */
	dt = dz; 
    } else {
	sf_putfloat(out,"d1",dt);
    }
    if (!sf_getfloat ("t0",&t0)) {
	/* Time origin (default is 0) */
	t0 = 0.; 
    } 
    sf_putfloat(out,"o1",t0);

    if (!sf_getbool ("slow",&slow)) slow = false;
    /* y: slowness, n: velocity */
    if (!sf_getfloat ("eps",&eps)) eps = 0.01;
    /* smoothness parameter */

    str = sf_stretch_init (nt, t0, dt, nz, eps, false);

    time = sf_floatalloc (nt);
    depth = sf_floatalloc (nz);
    vel = sf_floatalloc (nz);

    for (ix = 0; ix < nx; ix++) {
	sf_floatread (vel,nz,velocity);

	for (iz = 0; iz < nz; iz++) {
	    if (iz == 0) {
		z =  slow? z0*vel[0]/dz: z0/(dz*vel[0]);
	    } else {
		z += slow? vel[iz-1]: 1./vel[iz-1];
	    }
	    depth[iz] = 2.*dz*z;
	}

	sf_stretch_define (str, depth);

	sf_floatread (depth,nz,in);
	sf_stretch_apply (str, depth, time);
	sf_floatwrite (time,nt,out);
    }


    exit (0);
}

/* 	$Id: Mdepth2time.c 7107 2011-04-10 02:04:14Z ivlad $	 */
