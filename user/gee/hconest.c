/* Helical convolution, adjoint is the filter */
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

#include "hconest.h"

static float *x;
static sf_filter aa;

void hconest_init(float *x_in, sf_filter aa_in)
/*< initialize >*/
{
    x = x_in;
    aa = aa_in;
}

void hconest_lop(bool adj, bool add, int na, int ny, float *a, float *y)
/*< linear operator >*/
{
    int  ia, ix, iy;
    
    if (na != aa->nh) sf_error("%s: Wrong data dimensions",__FILE__);

    sf_adjnull(adj, add, na, ny, a, y);

    for (ia = 0; ia < na; ia++) {
	for (iy = aa->lag[ia]; iy < ny; iy++) {  
	    if(aa->mis[iy]) continue;
  
            ix = iy - aa->lag[ia];

	    if( adj) a[ia] -=  y[iy] * x[ix];
	    else     y[iy] -=  a[ia] * x[ix];
	}
    }
}

/* 	$Id: hconest.c 7107 2011-04-10 02:04:14Z ivlad $	 */
