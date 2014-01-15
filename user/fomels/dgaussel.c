/* Gaussian elimination (double precision) */
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
#include <float.h>

#include <rsf.h>

#include "dgaussel.h"

static int n;
static double **d;

void dgaussel_init (int size)
/*< initialize >*/
{
    int i;

    n = size;
    d = (double**) sf_alloc(n,sizeof(double*));
    d[0] = (double*) sf_alloc(n*(n+1),sizeof(double));
    for (i=1; i < n; i++) {
	d[i] = d[0]+i*(n+1);
    }
}

void dgaussel_close (void)
/*< free allocated storage >*/
{
    free(*d);
    free(d);
}

void dgaussel_solve (double **a       /* matrix [size][size] */, 
		     const double *b  /* rhs [size] */, 
		     float *x         /* solution [size] */)
/*< solve a*x = b >*/
{
    double dmax, di;
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
	dmax = fabs(d[k][k]);
	mloc = k;
	for (i=k+1; i < n; i++) {
	    di = fabs(d[i][k]);
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

	if (fabs(d[k][k]) < DBL_EPSILON) sf_error("%s: Bad matrix %g",__FILE__,d[k][k]);


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

