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
static sf_complex **o;

void cbanded_init (int n_in    /* matrix size */, 
		   int band_in /* band width */)
/*< Initialize >*/
{
    int i;

    n = n_in; 
    band = band_in;
    o = (sf_complex**) sf_alloc(band,sizeof(sf_complex*));
    for (i = 0; i < band; i++) {
	o[i] = sf_complexalloc (n-1-i);
    }
    d = sf_floatalloc(n);
}

void cbanded_const_define (float diag             /* diagonal */, 
			   const sf_complex *offd /* lower off-diagonal */)
/*< set matrix coefficients (constant along diagonals) >*/
{
    int k, m, j;
    sf_complex ct;
    float rt;

    d[0] = diag;
    for (k = 0; k < band-1; k++) {
	for (m = k; m >= 0; m--) {
	    ct = offd[m];
	    for (j = m+1; j < k-1; j++) {
#ifdef SF_HAS_COMPLEX_H
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
#else
		ct = sf_csub(ct,sf_cmul(o[j][k-j],
					sf_crmul(conjf(o[j-m-1][k-j]),
						 d[k-j])));
#endif
	    }
#ifdef SF_HAS_COMPLEX_H
	    o[m][k-m] = ct/d[k-m];
#else
	    o[m][k-m] = sf_crmul(ct,1./d[k-m]);
#endif
	}
	rt = diag;
	for (m = 0; m <= k; m++) {
#ifdef SF_HAS_COMPLEX_H
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
#else
	    rt -= crealf(sf_cmul(o[m][k-m],conjf(o[m][k-m])))*d[k-m];
#endif
	}
	d[k+1] = rt;
    }
    for (k = band-1; k < n-1; k++) {
	for (m = band-1; m >= 0; m--) {
	    ct = offd[m];
	    for (j = m+1; j < band; j++) {
#ifdef SF_HAS_COMPLEX_H
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
#else
		ct = sf_csub(ct,sf_cmul(o[j][k-j],
					sf_crmul(conjf(o[j-m-1][k-j]),
						 d[k-j])));
#endif
	    }
#ifdef SF_HAS_COMPLEX_H
	    o[m][k-m] = ct/d[k-m];
#else
	    o[m][k-m] = sf_crmul(ct,1./d[k-m]);
#endif
	}
	rt = diag;
	for (m = 0; m < band; m++) {
#ifdef SF_HAS_COMPLEX_H
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
#else
	     rt -= crealf(sf_cmul(o[m][k-m],conjf(o[m][k-m])))*d[k-m];
#endif
	}
	d[k+1] = rt;
    }
}

void cbanded_define (const float *diag    /* diagonal */, 
		     sf_complex **offd /* lower subdiagonals */)
/*< set matrix coefficients >*/
{
    int k, m, j;
    sf_complex ct;
    float rt;

    d[0] = diag[0];
    for (k = 0; k < band-1; k++) {
	for (m = k; m >= 0; m--) {
	    ct = offd[m][k-m];
	    for (j = m+1; j < k-1; j++) {
#ifdef SF_HAS_COMPLEX_H
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
#else
		ct = sf_csub(ct,sf_cmul(o[j][k-j],
					sf_crmul(conjf(o[j-m-1][k-j]),
						 d[k-j])));
#endif
	    }
#ifdef SF_HAS_COMPLEX_H
	    o[m][k-m] = ct/d[k-m];
#else
	    o[m][k-m] = sf_crmul(ct,1./d[k-m]);
#endif
	}
	rt = diag[k+1];
	for (m = 0; m <= k; m++) {
#ifdef SF_HAS_COMPLEX_H
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
#else
	    rt -= crealf(sf_cmul(o[m][k-m],conjf(o[m][k-m])))*d[k-m];
#endif
	}
	d[k+1] = rt;
    }
    for (k = band-1; k < n-1; k++) {
	for (m = band-1; m >= 0; m--) {
	    ct = offd[m][k-m];
	    for (j = m+1; j < band; j++) {
#ifdef SF_HAS_COMPLEX_H
		ct -= o[j][k-j]*conjf(o[j-m-1][k-j])*d[k-j];
#else
		ct = sf_csub(ct,sf_cmul(o[j][k-j],
					sf_crmul(conjf(o[j-m-1][k-j]),
						 d[k-j])));
#endif
	    }
#ifdef SF_HAS_COMPLEX_H
	    o[m][k-m] = ct/d[k-m];
#else
	    o[m][k-m] = sf_crmul(ct,d[k-m]);
#endif
	}
	rt = diag[k+1];
	for (m = 0; m < band; m++) {
#ifdef SF_HAS_COMPLEX_H
	    rt -= crealf(o[m][k-m]*conjf(o[m][k-m]))*d[k-m];
#else
	    rt -= crealf(sf_cmul(o[m][k-m],conjf(o[m][k-m])))*d[k-m];
#endif
	}
	d[k+1] = rt;
    }
}

void cbanded_solve (sf_complex *b)
/*< multiply by inverse (in place) >*/
{
    int k, m;
    sf_complex t;

    for (k = 1; k < band; k++) {
	t = b[k];
	for (m = 1; m <= k; m++) {
#ifdef SF_HAS_COMPLEX_H
	    t -= o[m-1][k-m] * b[k-m];
#else
	    t = sf_csub(t,sf_cmul(o[m-1][k-m],b[k-m]));
#endif
	}
	b[k] = t;
    }
    for (k = band; k < n; k++) {
	t = b[k];
	for (m = 1; m <= band; m++) {
#ifdef SF_HAS_COMPLEX_H
	    t -= o[m-1][k-m] * b[k-m];
#else
	    t = sf_csub(t,sf_cmul(o[m-1][k-m],b[k-m]));
#endif
	}
	b[k] = t;
    }
    for (k = n-1; k >= n - band; k--) {
#ifdef SF_HAS_COMPLEX_H
	t = b[k]/d[k];
#else
	t = sf_crmul(b[k],1./d[k]);
#endif
	for (m = 0; m < n - k - 1; m++) {
#ifdef SF_HAS_COMPLEX_H
	    t -= conjf(o[m][k]) * b[k+m+1];
#else
	    t = sf_csub(t,sf_cmul(conjf(o[m][k]),b[k+m+1]));
#endif
	}
	b[k] = t;
    }
    for (k = n - band - 1; k >= 0; k--) {
#ifdef SF_HAS_COMPLEX_H
	t = b[k]/d[k];
#else
	t = sf_crmul(b[k],1./d[k]);
#endif
	for (m = 0; m < band; m++) {
#ifdef SF_HAS_COMPLEX_H
	    t -= conjf(o[m][k]) * b[k+m+1];
#else
	    t = sf_csub(t,sf_cmul(conjf(o[m][k]),b[k+m+1]));
#endif
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

/* 	$Id: cbanded.c 7107 2011-04-10 02:04:14Z ivlad $	 */
