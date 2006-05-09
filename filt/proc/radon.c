/* Frequency-domain Radon transform */
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

#include "radon.h"

static float dp, p0;
static int nx, np;
static sf_complex *c0, *dc, czero;

void radon_init (int nx_in                            /* number of offsets */, 
		 int np_in , float dp_in, float p0_in /* slope axis */) 
/*< initialize >*/
{
    nx = nx_in; np = np_in;
    dp = dp_in; p0 = p0_in;
    c0 = sf_complexalloc(nx);
    dc = sf_complexalloc(nx);
    czero = sf_cmplx(0.,0.);
}

void radon_close () 
/*< free allocated storage >*/
{
    free (c0);
    free (dc);
}

void radon_set (float w   /* frequency */, 
		float* xx /* offset [nx] */)
/*< set the matrix >*/
{
    int ix;

    for (ix=0; ix < nx; ix++) {
      dc[ix] = sf_cmplx(cosf(w*dp*xx[ix]),sinf(w*dp*xx[ix]));
      c0[ix] = sf_cmplx(cosf(w*p0*xx[ix]),sinf(w*p0*xx[ix])); 
    }
}

void radon_toep (sf_complex *qq /* Toeplitz row */, 
		 float eps         /* regularization */)
/*< fill Toeplitz matrix for inversion >*/
{
    int ix, ip;
    sf_complex c, q;
    
    qq[0] = sf_cmplx(eps*eps,0.);
    for (ip=1; ip < np; ip++) {
	qq[ip] = czero;
    }

    for (ix=0; ix < nx; ix++) {
	c = conjf(dc[ix]);
	q = sf_cmplx(1.,0.);
	for (ip=0; ip < np-1; ip++) {
#ifdef SF_HAS_COMPLEX_H
	    qq[ip] += q;
	    q *= c;
#else
	    qq[ip] = sf_cadd(qq[ip],q);
	    q = sf_cmul(q,c);
#endif
	}
#ifdef SF_HAS_COMPLEX_H
	qq[np-1] += q;
#else
	qq[np-1] = sf_cadd(qq[np-1],q);
#endif
    }
}

void radon_lop (bool adj, bool add, int nm, int nd, 
		sf_complex *mm, sf_complex *dd)
/*< linear operator >*/
{
    int ix, ip;
    sf_complex c, d;

    if (nm != np || nd != nx) 
	sf_error("%s: mismatched data sizes",__FILE__);

    sf_cadjnull(adj, add, nm, nd, mm, dd);
    
    for (ix=0; ix < nx; ix++) {
	if (adj == true) {
	    c = dc[ix];
#ifdef SF_HAS_COMPLEX_H
	    d = dd[ix]*c0[ix];
#else
	    d = sf_cmul(dd[ix],c0[ix]);
#endif
	    for (ip=0; ip < np-1; ip++) {
#ifdef SF_HAS_COMPLEX_H
		mm[ip] += d;
		d *= c; 
#else
		mm[ip] = sf_cadd(mm[ip],d);
		d = sf_cmul(d,c); 
#endif
	    }
#ifdef SF_HAS_COMPLEX_H
	    mm[np-1] += d;
#else
	    mm[np-1] = sf_cadd(mm[np-1],d);
#endif
	} else {
	    c = conjf(dc[ix]);
	    d = mm[np-1];
	    for (ip=np-2; ip >= 0; ip--) {
#ifdef SF_HAS_COMPLEX_H
		d = d*c + mm[ip];
#else
		d = sf_cadd(sf_cmul(d,c),mm[ip]);
#endif
	    }
#ifdef SF_HAS_COMPLEX_H
	    dd[ix] += d*conjf(c0[ix]);
#else
	    dd[ix] = sf_cadd(dd[ix],sf_cmul(d,conjf(c0[ix])));
#endif
	}
    }
}

/* 	$Id$	 */

