/* Multiple arrivals by depth marching. */
/*
Copyright (C) 2008 University of Texas at Austin

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

#include "hdtrace.h"


int main(int argc, char* argv[])
{
    bool vel;
    int nz, nx, iz, np, npx, is, order, iorder, ix, k;
    float dz, dx, dp, x0, z0, p0;
    float **slow;
    float *slice[NS];
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

    /* default range from -1.0 to +1.0 (-90 to +90 from the vertical) */

    if (!sf_getint("np",&np)) np=201;
    /* number of horizontal slownesses */
    if (!sf_getfloat("dp",&dp)) dp=0.01;
    /* slowness increment (non dimensional) */
    if (!sf_getfloat("p0",&p0)) p0=-1.;
    /* starting horizontal slowness */

    for (is=0; is < NS; is++) {
        sf_putint(out[is],"n3",nz);
        sf_putfloat(out[is],"d3",dz);
        sf_putfloat(out[is],"o3",z0);

        sf_putint(out[is],"n2",nx);
        sf_putfloat(out[is],"d2",dx);
        sf_putfloat(out[is],"o2",x0);

        sf_putint(out[is],"n1",np);
        sf_putfloat(out[is],"d1",dp);
        sf_putfloat(out[is],"o1",p0);
    }

    npx = np*nx;

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
        slice[is] = sf_floatalloc(npx);
        for (k = 0; k < npx; k++) {
            slice[is][k] = 0.;
        }
    }

    hdtrace_init (order, iorder, nx, nz, np, dx, dz, dp, x0, z0, p0, slow, slice);

    /* loop for upgoing rays */
    for (iz = 0; iz < nz; iz++) {
	
	sf_warning("depth %d of %d", iz+1, nz);

	hdtrace_step (iz,-1,slow);

	for (is = 0; is < NS; is++) {
            sf_floatwrite (slice[is],npx,out[is]);
        }

    }

    exit (0);
}

/* 	$Id$	 */
