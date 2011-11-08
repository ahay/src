/* 2D Dynamic Ray Tracing */
/* NOTE: domain must be (0,z)*(0,x) with a central vertial ray at (0,0) */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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
    int dim, n[SF_MAX_DIM], i, j, s;
    float **der, **vel;
    float o[SF_MAX_DIM], d[SF_MAX_DIM], t0, shift, source, *t, dist;
    char key[6];
    sf_complex *P, *Q, **beam;
    sf_file in, deriv, out;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    /* read input dimension */
    dim = sf_filedims(in,n);

    if (dim > 2) sf_error("Only works for 2D.");

    for (i=0; i < dim; i++) {
	sprintf(key,"d%d",i+1);
	if (!sf_histfloat(in,key,d+i)) sf_error("No %s= in input.",key);
	sprintf(key,"o%d",i+1);
	if (!sf_histfloat(in,key,o+i)) o[i]=0.;
    }

    if (!sf_getfloat("t0",&t0)) t0=0.;
    /* time origin at source */
    
    if (!sf_getfloat("shift",&shift)) shift=1.;
    /* complex source shift */

    if (!sf_getfloat("source",&source)) source=o[1]+(n[1]-1)/2*d[1];
    /* source location */

    /* project source to grid point */
    s = (source-o[1])/d[1];
    
    /* read velocity model */
    vel = sf_floatalloc2(n[0],n[1]);
    sf_floatread(vel[0],n[0]*n[1],in);
    
    /* read derivative file */
    der = sf_floatalloc2(n[0],n[1]);
    if (NULL != sf_getstring("deriv")) {
	deriv = sf_input("deriv");
	sf_floatread(der[0],n[0]*n[1],deriv);
	sf_fileclose(deriv);
    } else {
	deriv = NULL;

	/* use finite difference for derivative if available */
	if (0 < s && s < n[1]-1) {
	    for (i=0; i < n[0]; i++)
		der[s][i] = (vel[s+1][i]-2*vel[s][i]+vel[s-1][i])/(d[1]*d[1]);
	} else {
	    for (i=0; i < n[0]; i++)
		der[s][i] = 0.;
	}
    }
    
    /* write output header */
    sf_settype(out,SF_COMPLEX);

    /* allocate memory for temporary data */
    t  = sf_floatalloc(n[0]);
    P  = sf_complexalloc(n[0]);
    Q  = sf_complexalloc(n[0]);
    beam = sf_complexalloc2(n[0],n[1]);

    /* complex source initial condition */
    t[0] = 0.;
    P[0] = sf_cmplx(0.,1./vel[s][0]);
    Q[0] = sf_cmplx(shift,0.);

    /* dynamic ray tracing along central ray (4th order Runge-Kutta) */
    for (i=1; i<n[0]; i++) {
	t[i] = t[i-1]+d[0]/2./vel[s][i-1]+d[0]/(vel[s][i-1]+vel[s][i]);
#ifdef SF_HAS_COMPLEX_H
	Q[i] = (P[i-1]-der[s][i-1]*Q[i-1]*d[0]/(vel[s][i-1]*vel[s][i-1]*4))
	    *vel[s][i-1]*d[0]+Q[i-1];
	P[i] = -((1.5*der[s][i-1]+0.5*der[s][i])*Q[i-1]+(der[s][i-1]+der[s][i])/2*vel[s][i-1]*d[0]/2*P[i-1])
	    *d[0]/(vel[s][i-1]*vel[s][i-1]*2)+P[i-1];
#else
	Q[i] = sf_cadd(sf_crmul(sf_cadd(P[i-1],sf_crmul(Q[i-1],-der[s][i-1]*d[0]/(vel[s][i-1]*vel[s][i-1]*4)))
				,vel[s][i-1]*d[0]),Q[i-1]);
	P[i] = sf_cadd(sf_crmul(sf_cadd(sf_crmul(Q[i-1],-((1.5*der[s][i-1]+0.5*der[s][i]))),
					sf_crmul(P[i-1],(der[s][i-1]+der[s][i])/2*vel[s][i-1]*d[0]/2)),d[0]/(vel[s][i-1]*vel[s][i-1]*2)),P[i-1]);
#endif
    }

    /* Gaussian beam */
    for (j=0; j<n[1]; j++) {
	dist = (j-s)*d[1];
	
	for (i=0; i<n[0]; i++) {
#ifdef SF_HAS_COMPLEX_H
	    beam[j][i] = t0+t[i]+0.5*dist*dist*P[i]/Q[i];
#else
	    beam[j][i] = sf_cadd(sf_cmplx(t0+t[i],0.),sf_crmul(sf_cdiv(P[i],Q[i]),0.5*dist*dist));
#endif
	}
    }

    sf_complexwrite(beam[0],n[0]*n[1],out);

    exit(0);
}
