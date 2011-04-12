/* 2-D Dynamic Ray Tracing */
/* Cheating: constant velocity approximation for central vertical ray */
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
    float **dir, **vel;
    float t0, s, z0, x0, dz, dx, *t;
    sf_file in, ve, cmplx;
    sf_complex *P, *Q, *cpxtbl;

    sf_init(argc,argv);
    in = sf_input("in");
    ve = sf_input("vel");
    cmplx = sf_output("out");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"o1",&z0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

    if (!sf_getfloat("t0",&t0)) t0=.0;
    if (!sf_getfloat("s",&s)) s=1.0;

    dir = sf_floatalloc2(nz,nx);
    vel = sf_floatalloc2(nz,nx);
    t = sf_floatalloc(nz);
    P = sf_complexalloc(nz);
    Q = sf_complexalloc(nz);
    cpxtbl = sf_complexalloc(nz*nx);

    sf_settype(cmplx,SF_COMPLEX);
    
    sf_floatread(dir[0],nz*nx,in);
    sf_floatread(vel[0],nz*nx,ve);
    sf_fileclose(ve);

    /* Complex source initial condition */
    t[0] = 0.;
    P[0] = sf_cmplx(0,1./vel[0][0]);
    Q[0] = sf_cmplx(s,0);

    /* Dynamic ray tracing along vertical ray */
    for (i=1; i<nz; i++) {
	t[i] = t[i-1]+dz/2./vel[0][i-1]+dz/(vel[0][i-1]+vel[0][i]);
#ifdef SF_HAS_COMPLEX_H
	Q[i] = (P[i-1]-dir[0][i-1]*Q[i-1]*dz/(vel[0][i-1]*vel[0][i-1]*4))*vel[0][i-1]*dz+Q[i-1];
	P[i] = -((1.5*dir[0][i-1]+0.5*dir[0][i])*Q[i-1]+(dir[0][i-1]+dir[0][i])/2*vel[0][i-1]*dz/2*P[i-1])*dz/(vel[0][i-1]*vel[0][i-1]*2)+P[i-1];
#else
	Q[i] = sf_cadd(sf_crmul(sf_cadd(P[i-1],sf_crmul(Q[i-1],-dir[0][i-1]*dz/(vel[0][i-1]*vel[0][i-1]*4))),vel[0][i-1]*dz),Q[i-1]);
	P[i] = sf_cadd(
	    sf_crmul(
		sf_cadd(sf_crmul(Q[i-1],-((1.5*dir[0][i-1]+0.5*dir[0][i]))),sf_crmul(P[i-1],(dir[0][i-1]+dir[0][i])/2*vel[0][i-1]*dz/2)),dz/(vel[0][i-1]*vel[0][i-1]*2)),P[i-1]);
#endif
    }

    /* Gaussian beam complex travel time */
    for (j=0; j<nx; j++) {
	for (i=0; i<nz; i++) {
#ifdef SF_HAS_COMPLEX_H
	    cpxtbl[i+j*nz] = t0+t[i]+0.5*(j*dx*j*dx)*P[i]/Q[i];	
#else
	    cpxtbl[i+j*nz] = sf_cadd(sf_cmplx(t0+t[i],0.),sf_crmul(sf_cdiv(P[i],Q[i]),0.5*(j*dx*j*dx)));	
#endif
	}
    }

    sf_complexwrite(cpxtbl,nz*nx,cmplx);

    exit(0);
}
