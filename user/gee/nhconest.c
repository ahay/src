/* Helical convolution for non-stationary filter estimation */
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

#include "nhconest.h"

#include "nhelix.h"
/*^*/

static float *x;
static nfilter aa;
static int nhmax;

void nhconest_init(float *x_in   /* data */, 
		   nfilter aa_in /* filter */, 
		   int nhmax_in  /* maximum filter size */)
/*< initialize >*/
{
    x = x_in;
    aa = aa_in;
    nhmax = nhmax_in;
}

void nhconest_lop(bool adj, bool add, int naa, int ny, float *a, float *y)
/*< linear operator >*/
{
    int ia, na, ix, iy, ip, *lag;

    if (naa != nhmax*aa->np) sf_error("%s: wrong dimensions",__FILE__);

    sf_adjnull(adj,add,naa,ny,a,y);

    for (iy=0; iy < ny; iy++) {
	if (aa->mis[iy]) continue;

        ip = aa->pch[iy]; 
	lag = aa->hlx[ip]->lag;
	na = aa->hlx[ip]->nh;

        for (ia=0; ia < na; ia++) {
	    ix = iy - lag[ia]; 
	    if (ix < 0) continue;
  
	    if (adj) {
		a[ia+nhmax*ip] += y[iy] * x[ix];
	    } else {
		y[iy] += a[ia+nhmax*ip] * x[ix];
	    }
        }
    }
}

