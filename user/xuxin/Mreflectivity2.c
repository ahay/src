/* 2-D reflectivity */
/*
  Copyright (C) 2011 KAUST
  
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

/* Adapted from Numerical Methods of Exploration Seismology by G Margrave */

#include <rsf.h>
#include <math.h>
#include "wave.h"

int main(int argc, char* argv[])
{
    bool dens;
    sf_file Fi,Fd,Fo;
    sf_axis ax,az;
    int nx,nz,ix,iz,n,order;
    float dx,dz,**v=NULL,**d=NULL,**r=NULL,**u=NULL,**ux=NULL,**uz=NULL,*t1=NULL,*t2=NULL;

    sf_init(argc,argv);
    Fi = sf_input("in");
    Fo = sf_output("out");
    if (NULL != sf_getstring("den")) {
	dens = true;
	Fd = sf_input("den");
    } else {
	dens = false;
	Fd = NULL;
    }
    
    if (!sf_getint("order",&order)) order=6;
    /* filter order */

    ax = sf_iaxa(Fi,1);
    az = sf_iaxa(Fi,2);
    nx = sf_n(ax); dx = sf_d(ax);
    nz = sf_n(az); dz = sf_d(az);
    if (dens) {
	if (!sf_histint(Fd,"n1",&n) || n != nx) sf_error("Need n1=%d in den",nx);
	if (!sf_histint(Fd,"n2",&n) || n != nz) sf_error("Need n2=%d in den",nz);
    }
    sf_setlabel(ax,"x"); sf_setunit(ax,"");
    sf_setlabel(az,"z"); sf_setunit(az,"");
    sf_oaxa(Fo,ax,1);
    sf_oaxa(Fo,az,2);

    /* velocity */
    v = sf_floatalloc2(nx,nz);
    sf_floatread(v[0],nx*nz,Fi);

    /* density */
    d = sf_floatalloc2(nx,nz);
    if (dens) {
	sf_floatread(d[0],nx*nz,Fd);
    } else {
	for (ix=0; ix < nx*nz; ix++)
	    d[0][ix] = 1.0f;
    }

    u = sf_floatalloc2(nx,nz);
    ux= sf_floatalloc2(nx,nz); /* grad(u) */
    uz= sf_floatalloc2(nx,nz);
    r = sf_floatalloc2(nx,nz); /* reflectivity */
    t1= sf_floatalloc(nz);
    t2= sf_floatalloc(nz);

    /* log(impedance) */
    for     (iz=0; iz<nz; iz++)
	for (ix=0; ix<nx; ix++)
	    u[iz][ix] = log(v[iz][ix]*d[iz][ix]);
/*
    diffx2d(u,ux,nx,nz,dx);
    diffz2d(u,uz,nx,nz,dz);
*/
    /* grad */
    sf_deriv_init(nx,order,0.);
    for (iz=0; iz < nz; iz++) {
	sf_deriv(u[iz],ux[iz]);
	for (ix=0; ix < nx; ix++)
	    ux[iz][ix] *= 1.0f/dx;
    }
    sf_deriv_free();
    sf_deriv_init(nz,order,0.);
    for (ix=0; ix < nx; ix++) {
	for (iz=0; iz<nz; iz++)
	    t1[iz] = u[iz][ix];
	sf_deriv(t1,t2);
	for (iz=0; iz<nz; iz++)
	    uz[iz][ix] = t2[iz] /dz;
    }
    sf_deriv_free();

    /* reflectivity */
    for     (iz=0; iz<nz; iz++)
	for (ix=0; ix<nx; ix++)
	    r[iz][ix] = 0.5*(ux[iz][ix]*ux[iz][ix]
                            +uz[iz][ix]*uz[iz][ix]);

    sf_floatwrite(r[0],nx*nz,Fo);

    exit(0);
}
