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
    int nz, nx, iz, np, npx, ix, order, iorder;
    float **slow, dz, dx, dp, x0, z0, p0, s;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    /* velocity or slowness */

    out = sf_output("out");
    /* escape time */

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
  
    if (!sf_histfloat(in,"o1",&z0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");

    if (!sf_getint("np",&np)) np=201;
    /* number of horizontal slownesses */
    if (!sf_getfloat("dp",&dp)) dp=0.01;
    /* slowness increment (non dimensional) */
    if (!sf_getfloat("p0",&p0)) p0=-1.;
    /* starting horizontal slowness */

    /* default range from -1.0 to +1.0 (-90 to +90 from the vertical) */

    sf_putint(out,"n4",nz);
    sf_putfloat(out,"d4",dz);
    sf_putfloat(out,"o4",z0);
    
    sf_putint(out,"n3",nx);
    sf_putfloat(out,"d3",dx);
    sf_putfloat(out,"o3",x0);
    
    sf_putint(out,"n2",np);
    sf_putfloat(out,"d2",dp);
    sf_putfloat(out,"o2",p0);

    sf_putint(out,"n1",NS+1);

    npx = np*nx;
    /* size of one depth slice */

    /* additional parameters */
    if(!sf_getbool("vel",&vel)) vel=true;
    /* y, input is velocity; n, slowness */

    if(!sf_getint("order",&order)) order=3;
    /* interpolation accuracy for velocity */
    if(!sf_getint("iorder",&iorder)) iorder=4;
    /* interpolation accuracy for grid */

    /* read velocity or slowness */
    slow  = sf_floatalloc2(nz,nx);
    sf_floatread(slow[0],nz*nx,in);

    /* convert to slowness */
    for(ix = 0; ix < nx; ix++){
	for (iz = 0; iz < nz; iz++) {
	    s = slow[ix][iz];
	    slow[ix][iz] = vel ? 1/s : s;
	}
    }

    hdtrace_init (order, iorder, nx, nz, np, dx, dz, dp, x0, z0, p0, slow);

    for (iz = 0; iz < nz; iz++) {
	sf_warning("depth %d of %d", iz+1, nz);
	hdtrace_step (iz,-1);
	hdtrace_write(out);
    }

    hdtrace_close();

    exit (0);
}

/* 	$Id$	 */
