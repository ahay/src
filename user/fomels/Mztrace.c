/* Multiple arrivals by depth marching. */
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

#include "ztrace.h"

int main(int argc, char* argv[])
{
    bool vel;
    int nz,nx, iz, na, nt, nax, ix, order, is, iorder, k;
    float **slow, *slice[NS];
    float dz,dx,da,x0,z0,a0;
    sf_file in, out[NS];

    sf_init(argc,argv);
    in = sf_input("in");
    out[0] = sf_output("out");
    out[1] = sf_output("place");
    out[2] = sf_output("depth");
    out[3] = sf_output("angle");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
  
    if (!sf_histfloat(in,"o1",&z0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_getint("nt",&nt)) nt=nx*nz;
    /* ray length bound */

    if (!sf_getint("na",&na)) na=362;
    /* number of angles */
    if (!sf_getfloat("da",&da)) da=0.5;
    /* angle increment (in degrees) */
    if (!sf_getfloat("a0",&a0)) a0=-90.;
    /* starting angle (in degrees) */

    for (is=0; is < NS; is++) {
	sf_putint(out[is],"n3",nz);
	sf_putfloat(out[is],"d3",dz);
	sf_putfloat(out[is],"o3",z0);

	sf_putint(out[is],"n2",nx);
	sf_putfloat(out[is],"d2",dx);
	sf_putfloat(out[is],"o2",x0);

	sf_putint(out[is],"n1",na);
	sf_putfloat(out[is],"d1",da);
	sf_putfloat(out[is],"o1",a0);
    }

    da *= (SF_PI/180.);
    a0 *= (SF_PI/180.);

    nax = na*nx;

    /* additional parameters */
    if(!sf_getbool("vel",&vel)) vel=true;
    /* y, input is velocity; n, slowness */
    if(!sf_getint("order",&order)) order=3;
    /* interpolation accuracy for velocity */
    if(!sf_getint("iorder",&iorder)) iorder=4;
    /* interpolation accuracy for grid */

    slow  = sf_floatalloc2(nz,nx);
    sf_floatread(slow[0],nz*nx,in);
    if (vel) { /* convert to slowness */
	for(ix = 0; ix < nx; ix++){
	    for (iz = 0; iz < nz; iz++) {
		slow[ix][iz] = 1./slow[ix][iz];
	    }
	}
    }
    
    for (is = 0; is < NS; is++) {
	slice[is] = sf_floatalloc(nax);
	for (k = 0; k < nax; k++) {
	    slice[is][k] = 0.;
	}
    }

    ztrace_init (order, iorder, nx, nz, na, nt,
                 dx, dz, da, x0, z0, a0, slow, slice);

    for (iz = 0; iz < nz; iz++) {
	sf_warning("depth %d of %d", iz+1, nz);

	ztrace_step (iz);

	for (is=0; is < NS; is++) {
	    sf_floatwrite (slice[is],nax,out[is]);
	}
    }

    exit (0);
}

/* 	$Id$	 */
