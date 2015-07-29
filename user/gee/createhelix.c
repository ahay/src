/* create a helix filter with a defined shape */
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

#ifdef __GNUC__
#ifndef alloca
#define alloca __builtin_alloca
#endif
#else /* not GNU C  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi) || defined(hpux) || defined(__hpux)
#include <alloca.h>
#endif
#endif

#include <rsf.h>

#include "createhelix.h"

sf_filter createhelix(int ndim    /* number of dimensions */, 
		      int* nd     /* data size [ndim] */, 
		      int* center /* filter center [ndim] */, 
		      int* gap    /* filter gap [ndim] */, 
		      int* na     /* filter size [ndim] */)
/*< allocate and output a helix filter >*/ 
{
    sf_filter aa;
    int ii[SF_MAX_DIM], na123, ia, nh, lag0a,lag0d, *lag, i;
    bool skip;

    for (na123 = 1, i=0; i < ndim; i++) na123 *= na[i];
    lag = (int*) alloca(na123*sizeof(int));

    /* index pointing to the "1.0" */
    lag0a = sf_cart2line (ndim, na, center); 

    nh=0;
    /* loop over linear index. */
    for (ia = 1+lag0a; ia < na123; ia++) { 
		sf_line2cart(ndim, na, ia, ii);
	
		skip = false;
		for (i=0; i < ndim; i++) {
			if (ii[i] < gap[i]) {
				skip = true;
				break;
			}
		}
		if (skip) continue;
	
		lag[nh] = sf_cart2line(ndim, nd, ii);
		nh++;                        /* got another live one */
    }
    /* center shift for nd cube */
    lag0d = sf_cart2line(ndim,  nd, center); 
    aa = sf_allocatehelix(nh); /* nh becomes size of filter */

    for (ia=0; ia < nh; ia++) {
		aa->lag[ia] = lag[ia] - lag0d; 
		aa->flt[ia] = 0.;
    }

    return aa;
}

/* 	$Id: createhelix.c 7267 2011-06-13 18:38:18Z saragiotis $	 */
