/* 2D Helmholtz solver by LU factorization. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
#include <umfpack.h>

#include "fdprep.h"

double maxvel(int nm, float *vel);

void fdprep(const double omega,
	    const float eps,
	    const int n1, const int n2,
	    const float d1, const float d2,
	    float **v,
	    const int npml,
	    const int pad1, const int pad2,
	    SuiteSparse_long n, SuiteSparse_long nz,
	    SuiteSparse_long *Ti, SuiteSparse_long *Tj,
	    double* Tx, double *Tz)
/*< 5-point finite-difference scheme >*/
{
    int i, j, index;
    double eta1, eta2, mvel, c1, c2;
    double *g1, *g2, **pad;
    double complex *s1, *s2, neib, cent;
    SuiteSparse_long count;
    
    /* prepare PML */
    eta1 = (double) npml*d1;
    eta2 = (double) npml*d2;
    
    mvel = maxvel(n1*n2,v[0]);
    c1 = -3.*log((double) eps)*mvel/(eta1*omega);
    c2 = -3.*log((double) eps)*mvel/(eta2*omega);

    g1 = (double*) sf_alloc(pad1,sizeof(double));
    for (i=0; i < pad1; i++) {
	if (i < npml) {
	    g1[i] = pow(c1*((npml-i)*d1/eta1),2.);
	} else if (i >= pad1-npml) {
	    g1[i] = pow(c1*((i-(pad1-npml-1))*d1/eta1),2.);
	} else {
	    g1[i] = 0.;
	}
    }

    s1 = (double complex*) sf_alloc(pad1,sizeof(double complex));
    for (i=0; i < pad1; i++) {
	s1[i] = 1./(1.+I*g1[i]);
    }

    g2 = (double*) sf_alloc(pad2,sizeof(double));
    for (j=0; j < pad2; j++) {
	if (j < npml) {
	    g2[j] = pow(c2*((npml-j)*d2/eta2),2.);
	} else if (j >= pad2-npml) {
	    g2[j] = pow(c2*((j-(pad2-npml-1))*d2/eta2),2.);
	} else {
	    g2[j] = 0.;
	}
    }

    s2 = (double complex*) sf_alloc(pad2,sizeof(double complex));
    for (j=0; j < pad2; j++) {
	s2[j] = 1./(1.+I*g2[j]);
    }

    /* extend model */
    pad = (double**) sf_alloc(pad2,sizeof(double*));
    pad[0] = (double*) sf_alloc(pad1*pad2,sizeof(double));
    for (j=1; j < pad2; j++) {
	pad[j] = pad[0]+j*pad1;
    }

    for (j=0; j < npml; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[0][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[0][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[0][n1-1];
	}
    }
    for (j=npml; j < pad2-npml; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[j-npml][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[j-npml][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[j-npml][n1-1];
	}
    }
    for (j=pad2-npml; j < pad2; j++) {
	for (i=0; i < npml; i++) {
	    pad[j][i] = v[n2-1][0];
	}
	for (i=npml; i < pad1-npml; i++) {
	    pad[j][i] = v[n2-1][i-npml];
	}
	for (i=pad1-npml; i < pad1; i++) {
	    pad[j][i] = v[n2-1][n1-1];
	}
    }

    /* assemble matrix in triplet form */
    count = 0;
    for (j=1; j < pad2-1; j++) {
	for (i=1; i < pad1-1; i++) {
	    index = (j-1)*(pad1-2)+(i-1);
	    
	    cent = 0.+I*0.;

	    /* left */
	    neib = conj((s1[i]/s2[j]+s1[i-1]/s2[j])/(2.*d1*d1));
	    cent += -neib;

	    if (i != 1) {
		Ti[count] = index;
		Tj[count] = index-1;
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* right */
	    neib = conj((s1[i]/s2[j]+s1[i+1]/s2[j])/(2.*d1*d1));
	    cent += -neib;

	    if (i != pad1-2) {
		Ti[count] = index;
		Tj[count] = index+1;
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* down */
	    neib = conj((s2[j]/s1[i]+s2[j-1]/s1[i])/(2.*d2*d2));
	    cent += -neib;

	    if (j != 1) {
		Ti[count] = index;
		Tj[count] = index-(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* up */
	    neib = conj((s2[j]/s1[i]+s2[j+1]/s1[i])/(2.*d2*d2));
	    cent += -neib;

	    if (j != pad2-2) {
		Ti[count] = index;
		Tj[count] = index+(pad1-2);
		Tx[count] = creal(neib);
		Tz[count] = cimag(neib);

		count++;
	    }

	    /* center */
	    cent += conj(pow(omega/pad[j][i],2.)/(s1[i]*s2[j]));
	    
	    Ti[count] = index;
	    Tj[count] = index;
	    Tx[count] = creal(cent);
	    Tz[count] = cimag(cent);
	    
	    count++;
	}
    }
}

void fdpad(const int npml,
	   const int pad1, const int pad2,
	   sf_complex **dat,
	   double *Bx, double *Bz)
/*< pad >*/
{
    int i, j;

    for (j=1; j < pad2-1; j++) {
	for (i=1; i < pad1-1; i++) {
	    if (i < npml || i >= pad1-npml || 
		j < npml || j >= pad2-npml) {
		Bx[(j-1)*(pad1-2)+(i-1)] = 0.;
		Bz[(j-1)*(pad1-2)+(i-1)] = 0.;
	    } else {
		Bx[(j-1)*(pad1-2)+(i-1)] = creal(dat[j-npml][i-npml]);
		Bz[(j-1)*(pad1-2)+(i-1)] = cimag(dat[j-npml][i-npml]);
	    }
	}
    }
}

void fdcut(const int npml,
	   const int pad1, const int pad2,
	   sf_complex **dat,
	   double *Xx, double *Xz)
/*< cut >*/
{
    int i, j;

    for (j=npml; j < pad2-npml; j++) {
	for (i=npml; i < pad1-npml; i++) {
	    dat[j-npml][i-npml] = sf_cmplx((float) Xx[(j-1)*(pad1-2)+(i-1)], 
					   (float) Xz[(j-1)*(pad1-2)+(i-1)]);
	}	    
    }
}

double maxvel(int nm, float *vel)
/* find maximum velocity */
{
    int i;
    float val;
    
    val = vel[0];
    for (i=1; i < nm; i++) {
	if (vel[i] > val) val = vel[i];
    }

    return ((double) val);
}
