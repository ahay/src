/* Complex Toeplitz solver */
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

#include "ctoeplitz.h"
/*^*/

static int n;
static sf_complex *a;

static sf_complex cdprod (int j, 
			  const sf_complex* a, const sf_complex* b) 
/* complex dot product */
{
    int i;
    sf_complex c;
    c = sf_cmplx(0.,0.);
    for (i=1; i <= j; i++) {
#ifdef SF_HAS_COMPLEX_H
	c += a[j-i]*conjf(b[i]);
#else
	c = sf_cadd(c,sf_cmul(a[j-i],conjf(b[i])));
#endif
    }
    return c;
}

void ctoeplitz_init (int n_in /* matrix size */)
/*< initialize >*/
{
    n = n_in;
    a = sf_complexalloc (n);
    a[0] = sf_cmplx(1.,0.);
}

void ctoeplitz_solve (const sf_complex *r /* top row of the matrix */, 
		      sf_complex *f       /* inverted in place */)
/*< apply the solver >*/
{    
    int i,j;
    sf_complex e,c,w, bot;
    float v;
    
    v=crealf(r[0]);
#ifdef SF_HAS_COMPLEX_H
    f[0] /= v;
#else
    f[0] = sf_crmul(f[0],1./v);
#endif
    
    for (j=1; j < n; j++) {
	e = cdprod(j,a,r);
#ifdef SF_HAS_COMPLEX_H
	c = -e/v;
#else
	c = sf_crmul(e,-1./v);
#endif

	v += crealf(c)*crealf(e) + cimagf(c)*cimagf(e);
       
	for (i=1; i <= j/2; i++) {
#ifdef SF_HAS_COMPLEX_H
	    bot  = a[j-i] + c*conjf(a[i]);
	    a[i] += c*conjf(a[j-i]);
#else
	    bot  = sf_cadd(a[j-i],sf_cmul(c,conjf(a[i])));
	    a[i] = sf_cadd(a[i],sf_cmul(c,conjf(a[j-i])));
#endif
	    a[j-i] = bot;
	}
	a[j] = c;
       
	w = cdprod(j,f,r);
#ifdef SF_HAS_COMPLEX_H
	c = (f[j]-w)/v;
#else
	c = sf_crmul(sf_csub(f[j],w),1./v);
#endif
       
	for (i=0; i < j; i++) {
#ifdef SF_HAS_COMPLEX_H
	    f[i] += c*conjf(a[j-i]);
#else
	    f[i] = sf_cadd(f[i],sf_cmul(c,conjf(a[j-i])));
#endif
	}
	f[j] = c;
    }
}

void ctoeplitz_close()
/*< free allocated storage >*/
{
    free (a);
}

/* 	$Id: ctoeplitz.c 7107 2011-04-10 02:04:14Z ivlad $	 */

