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

#include "helix_icai.h"

/*^*/

static sf_filter aa;
static int lg;

void helix_icai_init( sf_filter bb, int lag) 
/*<  Initialized with the filter. >*/
{
    aa = bb;
	lg = lag;
}

void  helix_icai_lop( bool adj, bool add, 
		     int nx, int ny, float* xx, float*yy) 
/*< linear operator >*/
{
    int ia, iy, ix;
    
	if(ny != nx) sf_error("%s: size problem: %d != %d",__FILE__,ny,nx);
	sf_adjnull(adj, add, nx, nx, xx, yy);

    for (ia = 0; ia < aa->nh; ia++) {
	
	for (iy = (aa->nh - lg) + aa->lag[ia]; iy <= ny-lg; iy++) {
	    
		if( aa->mis != NULL && aa->mis[iy]) continue;
	    ix = iy - aa->lag[ia] + lg - 1;
	    if(adj) {
		xx[ix] += yy[iy] * aa->flt[ia];
	    } else {
		yy[iy] += xx[ix] * aa->flt[ia];
	    }
	}
    }
}

/* 	$Id: helicon.c 2523 2007-02-02 16:45:29Z sfomel $	 */
