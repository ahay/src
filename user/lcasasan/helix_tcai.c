/* Helical transient convolution, adjoint is the input */
/*
  Copyright (C) 2010 Politecnico di Milano
  
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

#include "helix_tcai.h"


static sf_filter aa;

void helix_tcai_init( sf_filter aa_in)
/*< initialize >*/
{
    aa = aa_in;
}

void helix_tcai_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
/*< linear operator >*/
{
    int  ia, ix, iy;
    int na;
    
    na = aa->nh;
    //if (ny < nx+na-1) sf_error("%s: Wrong data dimensions ny=%d must be %d ",__FILE__,ny,nx+na-1);

    sf_adjnull(adj, add, nx, ny, x, y);

    for (ia = 0; ia < na; ia++) {
	for (iy = aa->lag[ia]; iy < ny; iy++) {
		if( aa->mis != NULL && aa->mis[iy]) continue;

            ix = iy - aa->lag[ia];

	    if( adj) x[ix] +=  y[iy] * aa->flt[ia];
	    else     y[iy] +=  x[ix] * aa->flt[ia];
	}
    }
}

/* 	$Id: hconest.c 2521 2007-02-02 00:25:42Z sfomel $	 */
