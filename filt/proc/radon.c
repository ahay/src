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
static float complex *c0, *dc, czero;

void radon_init (int nx_in                            /* number of offsets */, 
		 int np_in , float dp_in, float p0_in /* slope axis */) 
/*< initialize >*/
{
    nx = nx_in; np = np_in;
    dp = dp_in; p0 = p0_in;
    c0 = sf_complexalloc(nx);
    dc = sf_complexalloc(nx);
    czero = 0.;
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
      dc[ix] = cexpf(w*dp*xx[ix]*I);
      c0[ix] = cexpf(w*p0*xx[ix]*I); 
    }
}

void radon_toep (float complex *qq /* Toeplitz row */, 
		 float eps         /* regularization */)
/*< fill Toeplitz matrix for inversion >*/
{
    int ix, ip;
    float complex c, q;
    
    qq[0] = eps*eps;
    for (ip=1; ip < np; ip++) {
	qq[ip] = czero;
    }

    for (ix=0; ix < nx; ix++) {
	c = conjf(dc[ix]);
	q = 1.;
	for (ip=0; ip < np-1; ip++) {
	    qq[ip] += q;
	    q *= c;
	}
	qq[np-1] += q;
    }
}

void radon_lop (bool adj, bool add, int nm, int nd, 
		float complex *mm, float complex *dd)
/*< linear operator >*/
{
    int ix, ip;
    float complex c, d;

    if (nm != np || nd != nx) 
	sf_error("%s: mismatched data sizes",__FILE__);

    sf_cadjnull(adj, add, nm, nd, mm, dd);
    
    for (ix=0; ix < nx; ix++) {
	if (adj == true) {
	    c = dc[ix];
	    d = dd[ix]*c0[ix];
	    for (ip=0; ip < np-1; ip++) {
		mm[ip] += d;
		d *= c; 
	    }
	    mm[np-1] += d;
	} else {
	    c = conjf(dc[ix]);
	    d = mm[np-1];
	    for (ip=np-2; ip >= 0; ip--) {
		d = d*c + mm[ip];
	    }
	    dd[ix] += d*conjf(c0[ix]);
	}
    }
}

/* 	$Id$	 */

