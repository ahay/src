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

#include "ztrace2.h"

int main(int argc, char* argv[])
{
    bool vel;
    int nz,nx, iz, na, nax, ix, order, iorder;
    float **slow, dz,dx,da,x0,z0,a0;
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

    if (!sf_getint("na",&na)) na=362;
    /* number of angles */
    if (!sf_getfloat("da",&da)) da=0.5;
    /* angle increment (in degrees) */
    if (!sf_getfloat("a0",&a0)) a0=-90.;
    /* starting angle (in degrees) */

    /* default range from -90 to +90 from the vertical */

    sf_putint(out,"n4",nz);
    sf_putfloat(out,"d4",dz);
    sf_putfloat(out,"o4",z0);
    
    sf_putint(out,"n3",nx);
    sf_putfloat(out,"d3",dx);
    sf_putfloat(out,"o3",x0);
    
    sf_putint(out,"n2",na);
    sf_putfloat(out,"d2",da);
    sf_putfloat(out,"o2",a0);

    sf_putint(out,"n1",NS+1);

    /* convert degrees to radians */
    da *= (SF_PI/180.);
    a0 *= (SF_PI/180.);

    nax = na*nx;
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
    if (vel) {
	for(ix = 0; ix < nx; ix++){
	    for (iz = 0; iz < nz; iz++) {
		slow[ix][iz] = 1./slow[ix][iz];
	    }
	}
    }

    ztrace2_init (order, iorder, slow, 
		  nx, nz, na,
		  dx, dz, da, 
		  x0, z0, a0);

    for (iz = 0; iz < nz; iz++) {
	sf_warning("depth %d of %d", iz+1, nz);
	ztrace2_step (iz);
	ztrace2_write(out);
    }

    ztrace2_close();

    exit (0);
}

/* 	$Id$	 */
