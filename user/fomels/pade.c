/* Steering filter using Pade approximation */
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
#include <rsf.h>

#include "pade.h"

static float *w1, **w;
static int m;

void pade_init (int order)
/*< initialize >*/
{
    m = order;
    w = sf_floatalloc2(m,m);
    w1 = sf_floatalloc(m);
    sf_gaussel_init(m);
}    

void pade_close (void)
/*< free allocated storage >*/
{
    free(*w);
    free(w);
    free(w1);
    sf_gaussel_close();
}

void pade_zip (int n, const float *d, float *c)
/*< zip >*/
{
    int *k, i, j;

    k = sf_intalloc(n);
    for (i=0; i < n; i++) {
	k[i] = 0;
	c[i] = 0.;
    }
    k[0] = 1;

    for (i=0; i < n; i++) {
	for (j=i; j > 0; j--) {
	    k[j] += k[j-1];
	}
	for (j=0; j <= i; j++) {
	    c[j] += d[i]*k[j];
	}
    }
    free(k);
}

void pade_unzip (int n, float *d)
/*< unzip >*/
{
    float *c;
    int *k, i, j, s;

    k = sf_intalloc(n);
    c = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	k[i] = 0;
	c[i] = 0.;
    }
    k[0] = 1;

    for (i=0, s=1; i < n; i++, s=-s) {
	if (i > 0) {
	    for (j=i; j > 0; j--) {
		k[j] -= k[j-1];
	    }
	}
	for (j=0; j <= i; j++) {
	    c[j] += s*d[i]*k[j];
	}
    }

    for (i=0; i < n; i++) {
	d[i] = c[i];
    }

    free(k);
    free(c);
}

void pade_apply (int n, const float *c /* [n] */, int na, float * a, float *b)
/*< apply >*/
{
    int i, j;

    n -= m;
    if (na < n || m + 1 > n) sf_error("Wrong dimensions");

    /* form the matrix */
    for (i=0; i < m; i++) {
	for (j=0; j < m; j++) {
	    w[i][j] = c[n-i+j-1];
	}
	w1[i] = - c[n+i];
    }

    /* actually, a toeplitz matrix, but we are too lazy */
    sf_gaussel_solve (w,w1,b);

    /* convolution */
    for (i=0; i < n; i++) {
	a[i] = c[i];
    }
    for (i=0; i < m; i++) {
	for (j=i+1; j < na; j++) {
	    a[j] += b[i] * c[n-i+j-na-1];
	}
    }
}




