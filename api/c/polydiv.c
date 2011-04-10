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
#include <stdlib.h>

#include "polydiv.h"
#include "alloc.h"
#include "adjnull.h"

#include "helix.h"
/*^*/

static sf_filter aa;
static float* tt;

void sf_polydiv_init( int nd       /* data size */, 
		      sf_filter bb /* filter */) 
/*< initialize >*/
{
    aa = bb;
    tt = sf_floatalloc (nd);
}

void sf_polydiv_lop( bool adj, bool add, 
		     int nx, int ny, float* xx, float*yy) 
/*< linear operator >*/
{
    int ia, iy, ix;
    
    sf_adjnull( adj, add, nx, ny, xx, yy);
    
    for (ix=0; ix < nx; ix++) tt[ix] = 0.;
    
    if (adj) {
	for (ix = nx-1; ix >= 0; ix--) {  
	    tt[ix] = yy[ix];
	    for (ia = 0; ia < aa->nh; ia++) {
		iy = ix + aa->lag[ia];
		if( iy >= ny) continue;
		tt[ix] -= aa->flt[ia] * tt[iy];
	    }
	}
	for (ix=0; ix < nx; ix++) xx[ix] += tt[ix];
    } else {
	for (iy = 0; iy < ny; iy++) { 
	    tt[iy] = xx[iy];
	    for (ia = 0; ia < aa->nh; ia++) {
		ix = iy - aa->lag[ia];
		if( ix < 0) continue;
		tt[iy] -= aa->flt[ia] * tt[ix];
	    }
	}
	for (iy=0; iy < ny; iy++) yy[iy] += tt[iy];
    }
}

void sf_polydiv_close (void) 
/*< free allocated storage >*/
{
    free (tt);
}

/* 	$Id$	 */

