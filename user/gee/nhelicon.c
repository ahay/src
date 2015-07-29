/* Helical convolution with a non-stationary filter */
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

#include "nhelicon.h"

#include "nhelix.h"
/*^*/

static nfilter aa;

void nhelicon_init(nfilter aa_in)
/*< initialize with the filter >*/
{
    aa = aa_in;
}

void nhelicon_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int iy, ix, ia, ip, na, *lag;
    float *flt;

    sf_copy_lop(adj,add,nx,ny,xx,yy);

    for (iy=0; iy < ny; iy++) {    
        if (NULL != aa->mis && aa->mis[iy]) continue;
        ip = aa->pch[iy]; 
	lag = aa->hlx[ip]->lag; 
	flt = aa->hlx[ip]->flt;
	na = aa->hlx[ip]->nh;
        for (ia=0; ia < na; ia++) { 
            ix = iy - lag[ia]; 
	    if(ix < 0) continue;
            if (adj) {
		xx[ix] += yy[iy] * flt[ia];
	    } else {
		yy[iy] += xx[ix] * flt[ia];
            }
	}
    }
}
