/* Complex banded (Hermitian positive definite) matrix inversion */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
/*^*/

#include "cbanded.h"

static int n, band;
static float *d;
static float complex **o;

void cbanded_init (int n_in    /* matrix size */, 
		   int band_in /* band width */)
/*< Initialize >*/
{
    int i;

    n = n_in; 
    band = band_in;
    o = (float complex**) sf_alloc(band,sizeof(float complex*));
    for (i = 0; i < band; i++) {
	o[i] = sf_complexalloc (n-1-i);
    }
    d = sf_floatalloc(n);
}

void cbanded_const_define (float diag                /* diagonal */, 
			   const float complex *offd /* lower off-diagonal */)
/*< set matrix coefficients (constant along diagonals) >*/
{
    int k, m, j;
    float complex ct;
    float rt;

    d[0] = diag;
    for (k = 0; k < band-1; k++) {
	for (m = k; m >= 0; m--) {
	    ct = offd[m];
	    for (j = m+1; j < k-1; j++) {
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
	    }
	    o[m][k-m] = ct/d[k-m];
	}
	rt = diag;
	for (m = 0; m <= k; m++) {
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
	}
	d[k+1] = rt;
    }
    for (k = band-1; k < n-1; k++) {
	for (m = band-1; m >= 0; m--) {
	    ct = offd[m];
	    for (j = m+1; j < band; j++) {
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
	    }
	    o[m][k-m] = ct/d[k-m];
	}
	rt = diag;
	for (m = 0; m < band; m++) {
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
	}
	d[k+1] = rt;
    }
}

void cbanded_define (const float *diag    /* diagonal */, 
		     float complex **offd /* lower subdiagonals */)
/*< set matrix coefficients >*/
{
    int k, m, j;
    float complex ct;
    float rt;

    d[0] = diag[0];
    for (k = 0; k < band-1; k++) {
	for (m = k; m >= 0; m--) {
	    ct = offd[m][k-m];
	    for (j = m+1; j < k-1; j++) {
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
	    }
	    o[m][k-m] = ct/d[k-m];
	}
	rt = diag[k+1];
	for (m = 0; m <= k; m++) {
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
	}
	d[k+1] = rt;
    }
    for (k = band-1; k < n-1; k++) {
	for (m = band-1; m >= 0; m--) {
	    ct = offd[m][k-m];
	    for (j = m+1; j < band; j++) {
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
	    }
	    o[m][k-m] = ct/d[k-m];
	}
	rt = diag[k+1];
	for (m = 0; m < band; m++) {
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
	}
	d[k+1] = rt;
    }
}

void cbanded_solve (float complex *b)
/*< multiply by inverse (in place) >*/
{
    int k, m;
    float complex t;

    for (k = 1; k < band; k++) {
	t = b[k];
	for (m = 1; m <= k; m++) {
	    t -= o[m-1][k-m] * b[k-m];
	}
	b[k] = t;
    }
    for (k = band; k < n; k++) {
	t = b[k];
	for (m = 1; m <= band; m++) {
	    t -= o[m-1][k-m] * b[k-m];
	}
	b[k] = t;
    }
    for (k = n-1; k >= n - band; k--) {
	t = b[k]/d[k];
	for (m = 0; m < n - k - 1; m++) {
	    t -= conjf(o[m][k]) * b[k+m+1];
	}
	b[k] = t;
    }
    for (k = n - band - 1; k >= 0; k--) {
	t = b[k]/d[k];
	for (m = 0; m < band; m++) {
	    t -= conjf(o[m][k]) * b[k+m+1];
	}
	b[k] = t;
    }
}

void cbanded_close (void)
/*< free allocated storage >*/
{
    int i;

    for (i = 0; i < band; i++) {
	free(o[i]);
    }
    free (d);
    free (o);
}

/* 	$Id$	 */
