/* Helical convolution (adjoint is the input). */
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
#include <stdlib.h>
#include <rsf.h>

#include "helix_icaf.h"
/*^*/

static float *x;
static sf_filter aa;
static int lg;

void helix_icaf_init( float *x_in, sf_filter bb, int lag)
/*<  Initialized with the filter. >*/
{
    x = x_in;
	aa = bb;
	lg = lag;
}

void  helix_icaf_lop( bool adj, bool add, 
		     int na, int ny, float* a, float *y)
/*< linear operator >*/
{
    int ia, iy, ix;
    
    if (na != aa->nh) sf_error("%s: Wrong data dimensions",__FILE__);
    sf_adjnull(adj, add, na, ny, a, y);

    for (ia = 0; ia < aa->nh; ia++) {
	
	for (iy = aa->nh + aa->lag[ia]; iy <= ny; iy++) {
	    
		if( aa->mis != NULL && aa->mis[iy]) continue;
	    ix = iy - aa->lag[ia] - 1;
	    if(adj) {
	    a[ia] += y[iy-lg]  *  x[ix];
	    } else {
		y[iy-lg] += a[ia]  *  x[ix];
	    }
	}
    }
}

/* 	$Id: helicon.c 2523 2007-02-02 16:45:29Z sfomel $	 */
