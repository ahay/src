/* Complex matrix multiplication operator. */
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

#include "cmatmult.h"  
#include "adjnull.h"

#include "_bool.h"
#include "komplex.h"
/*^*/

static sf_complex **bb;

void sf_cmatmult_init(sf_complex **bb_in)
/*< initialize matrix >*/
{
    bb = bb_in;
}

void sf_cmatmult_lop (bool adj, bool add, int nx, int ny, 
		      sf_complex *x, sf_complex *y)
/*< operator >*/
{
    int ix, iy;
    sf_cadjnull (adj,add,nx,ny,x,y);

    for (ix=0; ix < nx; ix++) {
	for (iy=0; iy < ny; iy++) {
#ifdef SF_HAS_COMPLEX_H
	    if (adj) {
		x[ix] +=conjf(bb[iy][ix])*y[iy];
	    } else {
		y[iy] += bb[iy][ix]*x[ix];
	    }
#else
	    if (adj) {
		x[ix] = sf_cadd(x[ix],sf_cmul(conjf(bb[iy][ix]),y[iy]));
	    } else {
		y[iy] = sf_cadd(y[iy],sf_cmul(bb[iy][ix],x[ix]));
	    }
#endif
	}
    }
}
