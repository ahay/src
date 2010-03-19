/* 2-D Dynamic Ray Tracing */
/* Cheating: constant velocity appximation for central vertical ray */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    int nz, nx, i, j;
    float **vel;
    float t0, v0, s, z0, x0, dz, dx;
    sf_file inp, cmplx;
    sf_complex *P, *Q, *cpxtbl;

    sf_init(argc,argv);
    inp = sf_input("in");
    cmplx = sf_output("out");

    if (!sf_histint(inp,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(inp,"o1",&z0)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o2",&x0)) sf_error("No o2= in input");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in input");

    if (!sf_getfloat("v0",&v0)) v0=1.0;
    if (!sf_getfloat("t0",&t0)) t0=.0;
    if (!sf_getfloat("s",&s)) s=1.0;

    vel = sf_floatalloc2(nz,nx);
    P = sf_complexalloc(nz);
    Q = sf_complexalloc(nz);
    cpxtbl = sf_complexalloc(nz*nx);

    sf_settype(cmplx,SF_COMPLEX);
    
    sf_floatread(vel[0],nz*nx,inp);

    /* Complex source initial condition */
    P[0] = sf_cmplx(0,1/v0);
    Q[0] = sf_cmplx(s,0);

    /* Dynamic ray tracing along vertical ray */
    for (i=1; i<nz; i++) {
	Q[i] = (P[i-1]-vel[(nx-1)/2][i-1]*Q[i-1]*dz/(v0*v0*4))*v0*dz+Q[i-1];
	P[i] = -((1.5*vel[(nx-1)/2][i-1]+0.5*vel[(nx-1)/2][i])*Q[i-1]+(vel[(nx-1)/2][i-1]+vel[(nx-1)/2][i])/2*v0*dz/2*P[i-1])*dz/(v0*v0*2)+P[i-1];
    }

    /* Gaussian beam complex travel time */
    for (j=0; j<nx; j++) {
	for (i=0; i<nz; i++) {

	    cpxtbl[i+j*nz] = t0+(i*dz)/v0+0.5*((j-(nx-1)/2)*dx)*((j-(nx-1)/2)*dx)*P[i]/Q[i];

	}
    }

    sf_complexwrite(cpxtbl,nz*nx,cmplx);

    exit(0);
}
