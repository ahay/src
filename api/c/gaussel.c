/* Gaussian elimination */
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

#include "gaussel.h"

#include "_defs.h"
#include "alloc.h"
#include "error.h"
#include "c99.h"

static int n;
static float **d;

void sf_gaussel_init (int size)
/*< initialize >*/
{
    n = size;
    d = sf_floatalloc2(n+1,n);
}

void sf_gaussel_close (void)
/*< free allocated storage >*/
{
    free(*d);
    free(d);
}

void sf_gaussel_solve (float **a       /* matrix [size][size] */, 
		       const float *b  /* rhs [size] */, 
		       float *x        /* solution [size] */)
/*< solve a*x = b >*/
{
    float dmax, di;
    int mloc, k, i, j;

    for (i=0; i < n; i++) {
	for (j=0; j < n; j++) {
	    d[i][j] = a[i][j];
	}
    }
    for (i=0; i < n; i++) {
	d[i][n] = b[i];
    }

    for (k=0; k < n; k++) {
	/* pivoting */
	dmax = fabsf(d[k][k]);
	mloc = k;
	for (i=k+1; i < n; i++) {
	    di = fabsf(d[i][k]);
	    if (di > dmax) {
		dmax = di;
		mloc = i;
	    }
	}

	if (k != mloc) {
	    for (i=k; i <= n; i++) {
		di = d[k][i];
		d[k][i] = d[mloc][i];
		d[mloc][i] = di;
	    }
	}

	if (fabsf(d[k][k]) < SF_EPS) sf_error("%s: Bad matrix",__FILE__);


	/* triangularization phase */
	for(j=k+1; j<n; j++) {
	    di = d[j][k]/d[k][k];
	    for(i=k; i<=n; i++) {
		d[j][i] -= di*d[k][i];
	    }
	}
    }

    /* back substitution phase */
    for (k=n-1; k >=0; k--) {
	di = d[k][n];
	for (i=k+1; i < n; i++) {
	    di -= d[k][i]*x[i];
	}
	x[k] = di/d[k][k];
    }
}

