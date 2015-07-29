/* Put a helical filter in a box. */
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

#include "boxfilter.h"

void box (int dim            /* number of dimaneions */, 
	  const int *nd      /* data size [dim] */, 
	  const int *center  /* filter center [dim] */, 
	  const int *na      /* filter size [dim] */, 
	  const sf_filter aa /* input filter */, 
	  int nc             /* box size */, 
	  float* cube        /* output box [nc] */)
/*< box it >*/ 
{
    int ii[SF_MAX_DIM];
    int j, lag0a, lag0d, id, ia;

    for (ia=0; ia < nc; ia++) {
	cube[ia] = 0.;
    }
    lag0a = sf_cart2line(dim, na, center);  /* 1.0 in na. */
    cube[lag0a] = 1.;                       /* place it. */
    lag0d = sf_cart2line(dim, nd, center);  /* 1.0 in nd. */
    for (j=0; j < aa->nh; j++) { /* inspect the entire helix */
	id = aa->lag[j] + lag0d;
	sf_line2cart(dim, nd, id, ii);	/* ii = cartesian  */
	ia = sf_cart2line(dim, na, ii);	/* ia = linear in aa */
	cube[ia] = aa->flt[j];		/* copy the filter */
    }
}

/* 	$Id: boxfilter.c 7107 2011-04-10 02:04:14Z ivlad $	 */

