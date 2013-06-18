/* Source for the method of manufactured solution */
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
#include <math.h>



int main(int argc, char* argv[])
{
    /*grid parameters*/
    int nx, nz, nt;
    float dx, dz, dt;
    int ix, iz, it;
    
    sf_axis at, ax, az;

    float **vel;          /*input velocity*/
    float **src, **slt; /*output MMS source and solution*/
    
    /*parameters*/
    float slx, slz;
    float alpha, beta, beta2;
    float dist2, xx, zz, tt, tmp;
    
    sf_file Fvel, Fsrc, Fslt;
    sf_init(argc, argv);
    Fvel = sf_input("in");
    Fsrc = sf_output("out");
    Fslt = sf_output("mslt"); /*the manufactured solution*/
    
    if (SF_FLOAT != sf_gettype(Fvel)) sf_error("Need float input");
    ax = sf_iaxa(Fvel, 2); nx = sf_n(ax); dx = sf_d(ax); 
    az = sf_iaxa(Fvel, 1); nz = sf_n(az); dz = sf_d(az);

    if (!sf_getint("nt",&nt)) sf_error("Need nt");
    /*number of time step*/
    if (!sf_getfloat("dt", &dt)) sf_error("Need dt");
    /*time step*/
    if (!sf_getfloat("slx", &slx)) slx = nx*dx*0.5;
    /*center of source location: x*/
    if (!sf_getfloat("slz", &slz)) slz = nz*dz*0.5;
    /*center of source location: z*/
    if (!sf_getfloat("alpha", &alpha)) alpha = 1.0e-2;
    /*source parameter*/
    if (!sf_getfloat("beta", &beta)) beta = 1.0;
    /*source parameter*/

    at = sf_maxa(nt, 0.0, dt);
        
    vel = sf_floatalloc2(nz, nx);
    src = sf_floatalloc2(nz, nx);
    slt = sf_floatalloc2(nz, nx);
    
    sf_floatread(vel[0], nx*nz, Fvel);

    /*Set output axis*/
    sf_oaxa(Fsrc, az, 1);
    sf_oaxa(Fsrc, ax, 2);
    sf_oaxa(Fsrc, at, 3);
    
    sf_oaxa(Fslt, az, 1);
    sf_oaxa(Fslt, ax, 2);
    sf_oaxa(Fslt, at, 3);

    /*Manufactured Solution Source*/
    beta2 = beta*beta;
    for (it=0; it<nt; it++) {
	//tt = cosf(it*dt*beta)/beta;
	tt = expf(-1*(it*dt-0.15)*beta);
	for (ix=0; ix<nx; ix++) {
	    xx = ix*dx;
	    for (iz=0; iz<nz; iz++) {
		zz = iz*dz;
		dist2 = (xx-slx)*(xx-slx) + (zz-slz)*(zz-slz);
		tmp = expf(-1*alpha*dist2)*tt;
		//src[ix][iz] = tmp*(beta2+4*alpha*vel[ix][iz]*vel[ix][iz]*(alpha*dist2-1)); 
		src[ix][iz] = tmp*(-1*beta2+4*alpha*vel[ix][iz]*vel[ix][iz]*(alpha*dist2-1)); 

	    }
	    sf_floatwrite(src[ix], nz, Fsrc);
	}
    }
    
    
    /*solution*/
    for (it=0; it<nt; it++) {
	//tt = sinf(it*dt*beta);
	tt = expf(-1*(it*dt-0.15)*beta);
	for (ix=0; ix<nx; ix++) {
	    xx = ix*dx;
	    for (iz=0; iz<nz; iz++) {
		zz = iz*dz;
		dist2 = (xx-slx)*(xx-slx) + (zz-slz)*(zz-slz);
		slt[ix][iz] = expf(-1*alpha*dist2)*tt;
	    }
	    sf_floatwrite(slt[ix], nz, Fslt);
	}
    }
    
}
					  
