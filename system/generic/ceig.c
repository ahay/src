/* Eigenvalues of a complex matrix. */
/*
  Copyright (C) 2016 The University of Texas at Austin

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

#include "jacobi2.h"

static bool verb;
static int n, n2;
static float *rwork;
static sf_complex *b, *work;

void ceig_init(bool verb1, int n1)
/*< allocate storage >*/
{
    verb = verb1;
    n = n1;
    n2 = 2*n;

    b = sf_complexalloc(n*n);
    work = sf_complexalloc(n2);
    rwork = sf_floatalloc(n2);

    jacobi2_init(n,verb);
}

void ceig_close(void)
/*< free allocated storage >*/
{
    free(b);
    free(work);
    free(rwork);

    jacobi2_close();
}


void ceig(int niter      /* number of iterations */, 
	  float tol      /* tolerance */, 
	  int m          /* effective matrix size */,
	  sf_complex** a /* [n][n] matrix */, 
	  sf_complex *e  /* [n] eigenvalues */)
/*< find eigenvalues >*/
{
    int iter, j, k, info;
    float s2,s0=1.;

    if (niter > 0) { /* Jacobi iterations */
	for (iter=0; iter < niter; iter++) {
	    s2 = 0.;
	    for (j=0; j < m; j++) {
		for (k=0; k < m; k++) {
		    s2 += jacobi2(a,m,j,k);
		}
	    }
	    if (verb) sf_warning("iter=%d s2=%g",iter+1,s2);
	    if (0==iter) {
		s0 = s2;
	    } else {
		if (s2 <= s0*tol) break;
	    }
	}
	
	for (j=0; j < m; j++) {
	    e[j]=a[j][j];
	}
    } else {
	for (j=0; j < m; j++) {
	    for (k=0; k < m; k++) {
		b[k+j*m] = a[j][k];
	    }
	}
#ifdef SF_HAS_LAPACK
	cgeev_( "N", "N", &m, b, &m, e, work, &n2, work, &n2, work, &n2, rwork, &info );
	if (info) sf_error("cgeev_ failed");
#else
	sf_error("No LAPACK");
#endif
    }
    for (j=m; j < n-1; j++) {
	e[j]=sf_cmplx(0.,0.);
    }
}
	
