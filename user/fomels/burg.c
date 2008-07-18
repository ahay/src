/* 1-D PEF estimation with complex data using Burg's algorithm */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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
/* The algorithm is described in J.F. Clarbout, FGDP, p. 133-... */

#include <rsf.h>

#include "burg.h"

static  int n, nf, nc;
float **f, **b;

void burg_init (int n_in  /* trace length */, 
		 int nc_in /* number of traces */, 
		 int nf_in /* filter length */)
/*< initialize >*/
{
    n = n_in; 
    nf = nf_in; 
    nc = nc_in;

    f = sf_floatalloc2(n,nc);
    b = sf_floatalloc2(n,nc);
}

void burg_close(void)
/*< free allocated storage >*/
{
    free (*f);
    free (f);
    free (*b);
    free (b);
}

void burg_apply (float *x  /* [n*nc] input data */, 
		 float *a  /* [nf] output prediction-error filter */)
/*< estimate PEF >*/
{
    double den, num, cj;
    float fi, bi, ai;
    int j, ic, i;

    for (ic=0; ic < nc; ic++) {
	for (i=0; i < n; i++) {
	    b[ic][i] = f[ic][i] = x[ic*n+i];
	}
    }

    a[0] = 1.;
    for (j=1; j < nf; j++) {
	num = den = 0.;
	for (ic=0; ic < nc; ic++) {
	    for (i=j; i < n; i++) {
		fi = f[ic][i];
		bi = b[ic][i-j];

		num += fi*bi;
		den += fi*fi + bi*bi;
	    }
	}
	cj = 2.*num/den;

	for (ic=0; ic < nc; ic++) {
	    for (i=j; i < n; i++) {
		fi = f[ic][i];
		bi = b[ic][i-j];

		f[ic][i] -= cj*bi;
		b[ic][i-j] -= cj*fi;
	    }
	}
	for (i=1; i <= j/2; i++) {
	    ai = a[j-i]-cj*a[i];
	    a[i] -= cj*a[j-i];
	    a[j-i] = ai;
	}
	a[j] = -cj;
    }
}
