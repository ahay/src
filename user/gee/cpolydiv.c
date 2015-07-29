/* Helical inverse convolution (polynomial division). */
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

#include "chelix.h"
#include "cpolydiv.h"

static cfilter aa;
static sf_complex* tt;

void cpolydiv_init( int nd       /* data size */, 
		   cfilter bb /* filter */) 
/*< initialize >*/
{
    aa = bb;
    tt = sf_complexalloc (nd);
}

void cpolydiv_lop( bool adj, bool add, 
		   int nx, int ny, sf_complex* xx, sf_complex* yy) 
/*< linear operator >*/
{
    int ia, iy, ix;
    
    sf_cadjnull( adj, add, nx, ny, xx, yy);
    
    for (ix=0; ix < nx; ix++) tt[ix] = sf_cmplx(0.,0.);
    
    if (adj) {
	for (ix = nx-1; ix >= 0; ix--) {  
	    tt[ix] = yy[ix];
	    for (ia = 0; ia < aa->nh; ia++) {
		iy = ix + aa->lag[ia];
		if( iy >= ny) continue;
#ifdef SF_HAS_COMPLEX_H
		tt[ix] -= conjf(aa->flt[ia]) * tt[iy];
#else
		tt[ix] = sf_cadd(tt[ix],
				 sf_cmul(conjf(aa->flt[ia]),sf_cneg(tt[iy])));
#endif
	    }
	}
	for (ix=0; ix < nx; ix++) { 
#ifdef SF_HAS_COMPLEX_H
	    xx[ix] += tt[ix];
#else
	    xx[ix] = sf_cadd(xx[ix],tt[ix]);
#endif
	}
    } else {
	for (iy = 0; iy < ny; iy++) { 
	    tt[iy] = xx[iy];
	    for (ia = 0; ia < aa->nh; ia++) {
		ix = iy - aa->lag[ia];
		if( ix < 0) continue;
#ifdef SF_HAS_COMPLEX_H
		tt[iy] -= aa->flt[ia] * tt[ix];
#else
		tt[iy] = sf_cadd(tt[iy],sf_cmul(aa->flt[ia],sf_cneg(tt[ix])));
#endif
	    }
	}
	for (iy=0; iy < ny; iy++) {
#ifdef SF_HAS_COMPLEX_H
	    yy[iy] += tt[iy];
#else
	    yy[iy] = sf_cadd(yy[iy],tt[iy]);
#endif
	}
    }
}

void cpolydiv_close (void) 
/*< free allocated storage >*/
{
    free (tt);
}



