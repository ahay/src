/* Interface for frequency-domain image interpolation */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

#include "freqmig0.h"

static int n1, n2, nf;
static float *omega;
static sf_complex **img, *x0;

void freqmig0_init(int nn1, int nn2              /* image dimensions */,
		   float of0, float df0, int nf0 /* frequency samples */,
		   sf_complex **iimg             /* data */)
/*< initialize >*/
{
    int iw;

    n1 = nn1; n2 = nn2; nf = nf0;
    img = iimg;

    /* frequency to omega */
    omega = sf_floatalloc(nf);

    for (iw=0; iw < nf; iw++) {
	omega[iw] = 2.*SF_PI*(of0+iw*df0);
    }
}

void freqmig0_close()
/*< free allocated storage >*/
{
    free(omega);
}

void freqmig0_set(sf_complex *xx0 /* current coefficients */)
/*< set-up >*/
{
    x0 = xx0;
}

float freqmig0_cost(sf_complex *coe /* coefficients */,
		    sf_complex *dr  /* right-hand side */)
/*< evaluate right-hand side and cost >*/
{
    int i, iw, nn;
    float cost=0.;
    sf_complex temp;

    for (i=0; i < n1*n2; i++) {
	nn = n1*n2+i;

	for (iw=0; iw < nf; iw++) {
	    temp = sf_cmplx(coe[i]*cosf(omega[iw]*coe[nn]),
			    coe[i]*sinf(omega[iw]*coe[nn]));
	    
	    /* right-hand side */
	    dr[n1*n2*iw+i] = img[iw][i]-temp;

	    /* least-squares cost */
	    cost += 0.5*crealf(conjf(temp-img[iw][i])*(temp-img[iw][i]));
	}
    }

    return cost;
}

void freqmig0_oper(bool adj, bool add, int nx, int nr, sf_complex *x, sf_complex *r)
/*< linear operator >*/
{
    int i, iw, nn;
    sf_complex temp;

    sf_cadjnull(adj,add,nx,nr,x,r);

    if (adj) {
	/* given dr solve dx */
	
	for (i=0; i < n1*n2; i++) {
	    nn = n1*n2+i;

	    for (iw=0; iw < nf; iw++) {
		/* A */
		temp = sf_cmplx(cosf(omega[iw]*x0[nn]),
				-sinf(omega[iw]*x0[nn]));
		x[i] += temp*r[n1*n2*iw+i];

		/* T */
		temp = sf_cmplx(-omega[iw]*x0[i]*sinf(omega[iw]*x0[nn]),
				-omega[iw]*x0[i]*cosf(omega[iw]*x0[nn]));
		x[nn] += temp*r[n1*n2*iw+i];
	    }
	}
    } else {
	/* given dx solve dr */

	for (i=0; i < n1*n2; i++) {
	    nn = n1*n2+i;

	    for (iw=0; iw < nf; iw++) {
		/* A */
		temp = sf_cmplx(cosf(omega[iw]*x0[nn]),
				sinf(omega[iw]*x0[nn]));
		r[n1*n2*iw+i] += temp*x[i];

		/* T */
		temp = sf_cmplx(-omega[iw]*x0[i]*sinf(omega[iw]*x0[nn]),
				omega[iw]*x0[i]*cosf(omega[iw]*x0[nn]));
		r[n1*n2*iw+i] += temp*x[nn];
	    }
	}
    }
}
