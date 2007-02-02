/* Compress a helical filter */
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

#include <math.h>

#include <rsf.h>
/*^*/

#include "compress.h"

sf_filter compress( sf_filter aa /* filter */, 
		    float eps    /* threshold */) 
/*< return a new filter with coefficients less than eps removed >*/
{
    sf_filter bb;
    bool* keep;
    int i, k;

    keep = sf_boolalloc (aa->nh);
    k = 0;
    for(i=0; i < aa->nh; i++) {
	if ( fabsf( aa->flt[i]) > eps) {
	    keep[i] = true;
	    k++;
	} else {
	    keep[i] = false;
	}
    }
    bb = sf_allocatehelix( k);
    k = 0;
    for(i=0; i < aa->nh; i++) {
	if (keep[i]) {
	    bb->flt[k] = aa->flt[i];
	    bb->lag[k] = aa->lag[i];
	    k++;
	}
    }
    sf_deallocatehelix( aa);
    free (keep);

    return bb;
}

/* 	$Id$	 */
