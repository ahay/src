/* Helical convolution. */
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

#include "helicon.h"
#include "copy.h"

#include "helix.h"
/*^*/

static sf_filter aa;

void sf_helicon_init( sf_filter bb) 
/*<  Initialized with the filter. >*/
{
    aa = bb;
}

void sf_helicon_lop( bool adj, bool add, 
		     int nx, int ny, float* xx, float*yy) 
/*< linear operator >*/
{
    int ia, iy, ix;
    
    sf_copy_lop(adj, add, nx, nx, xx, yy);

    if(adj) {
        for (ia = 0; ia < aa->nh; ia++) {
	    for (iy = aa->lag[ia]; iy < nx; iy++) {
	        if( aa->mis != NULL && aa->mis[iy]) continue;
	        ix = iy - aa->lag[ia];
	        xx[ix] += yy[iy] * aa->flt[ia];
	    }
	}
    } else {
        for (ia = 0; ia < aa->nh; ia++) {
	    for (iy = aa->lag[ia]; iy < nx; iy++) {
	        if( aa->mis != NULL && aa->mis[iy]) continue;
	        ix = iy - aa->lag[ia];
	        yy[iy] += xx[ix] * aa->flt[ia];
	    }
	}
    }
}

/* 	$Id$	 */
