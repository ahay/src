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

#include "cburg.h"

static  int n, nf, nc;
sf_complex **f, **b;

void cburg_init (int n_in  /* trace length */, 
		 int nc_in /* number of traces */, 
		 int nf_in /* filter length */)
/*< initialize >*/
{
    n = n_in; 
    nf = nf_in; 
    nc = nc_in;

    f = sf_complexalloc2(n,nc);
    b = sf_complexalloc2(n,nc);
}

void cburg_close(void)
/*< free allocated storage >*/
{
    free (*f);
    free (f);
    free (*b);
    free (b);
}

void cburg_apply (sf_complex *x  /* [n*nc] input data */, 
		  sf_complex *a  /* [nf] output prediction-error filter */)
/*< estimate PEF >*/
{
    double den;
    sf_complex fi, bi, ai, num, cj;
    int j, ic, i;

    for (ic=0; ic < nc; ic++) {
	for (i=0; i < n; i++) {
	    b[ic][i] = f[ic][i] = x[ic*n+i];
	}
    }

    a[0] = sf_cmplx(1.,0.);
    for (j=1; j < nf; j++) {
	num = sf_cmplx(0.,0.);
	den = 0.;
	for (ic=0; ic < nc; ic++) {
	    for (i=j; i < n; i++) {
		fi = f[ic][i];
		bi = b[ic][i-j];
#ifdef SF_HAS_COMPLEX_H
		num += fi*conj(bi);
		den += creal(fi*conj(fi) + bi*conj(bi));
#else
		num = sf_cadd(num,sf_cmul(fi,conjf(bi)));
		den += sf_crealf(sf_cadd(sf_cmul(fi,conjf(fi)),
					 sf_cmul(bi,conjf(bi))));
#endif
	    }
	}
#ifdef SF_HAS_COMPLEX_H
	cj = 2.*num/den;
#else
	cj = sf_crmul(num,2./den);
#endif
	for (ic=0; ic < nc; ic++) {
	    for (i=j; i < n; i++) {
		fi = f[ic][i];
		bi = b[ic][i-j];
#ifdef SF_HAS_COMPLEX_H
		f[ic][i] -= cj*bi;
		b[ic][i-j] -= conj(cj)*fi;
#else
		f[ic][i] = sf_cadd(f[ic][i],sf_cneg(sf_cmul(cj,bi)));
		b[ic][i-j] = sf_cadd(b[ic][i-j],
				     sf_cneg(sf_cmul(sf_conjf(cj),fi)));
#endif
	    }
	}
	for (i=1; i <= j/2; i++) {
#ifdef SF_HAS_COMPLEX_H
	    ai = a[j-i]-cj*conjf(a[i]);
	    a[i] -= cj*conjf(a[j-i]);
#else
	    ai = sf_cadd(a[j-i],sf_cneg(sf_cmul(cj,sf_conjf(a[i]))));
	    a[i] = sf_cadd(a[i],sf_cneg(sf_cmul(cj,sf_conjf(a[j-i]))));
#endif
	    a[j-i] = ai;
	}
#ifdef SF_HAS_COMPLEX_H
	a[j] = -cj;
#else
	a[j] = sf_cneg(cj);
#endif
    }
}
