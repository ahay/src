/* First-arrival traveltime table using analytical traveltimes */
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
#include <rsf.h>

#include "analytical.h"

int main(int argc, char* argv[])
{
    int nz,nx, iz,ix, order;
    bool vel;
    float x1,x2,z1,z2,v1,v2,x0,dx,z0,dz,fx,fz;
    float **slow, *time, g1[2],g2[2],p1[2],p2[2];
    sf_eno2 cvel;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    /* velocity or slowness */

    out = sf_output("out");
    /* traveltime */

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
  
    if (!sf_histfloat(in,"o1",&z0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");
    
    if(!sf_getbool("vel",&vel)) vel=true;
    /* y, input is velocity; n, slowness */

    if(!sf_getint("order",&order)) order=3;
    /* interpolation accuracy for velocity */

    if (!sf_getfloat("yshot",&x1)) x1=x0 + 0.5*(nx-1)*dx;
    if (!sf_getfloat("zshot",&z1)) z1=z0; 

    /* read velocity or slowness */
    slow  = sf_floatalloc2(nz,nx);
    sf_floatread(slow[0],nz*nx,in);    

    /* convert to slowness squared */
    for(ix = 0; ix < nx; ix++){
	for (iz = 0; iz < nz; iz++) {
	    v1 = slow[ix][iz];
	    v1 *= v1;
	    slow[ix][iz] = vel? 1/v1:v1;
	}
    }

    cvel = sf_eno2_init (order, nz, nx);
    sf_eno2_set (cvel, slow);

    time = sf_floatalloc(nz);

    fx = (x1-x0)/dx;
    ix = floorf(fx);
    fx -= ix;
    
    fz = (z1-z0)/dz;
    iz = floorf(fz);
    fz -= iz;

    sf_eno2_apply(cvel,iz,ix,fz,fx,&v1,g1,BOTH);
    g1[1] /= dx;
    g1[0] /= dz;
    
    for(ix = 0; ix < nx; ix++){
	x2 = x0+ix*dx;
	for (iz = 0; iz < nz; iz++) {
	    z2 = z0+iz*dz;

	    sf_eno2_apply(cvel,iz,ix,0.,0.,&v2,g2,BOTH);
	    g2[1] /= dx;
	    g2[0] /= dz;

	    time[iz] = analytical(x1,z1,v1,g1,p1,
				  x2,z2,v2,g2,p2);
	}
	sf_floatwrite(time,nz,out);
    }

    exit(0);
}
