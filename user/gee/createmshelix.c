/* Create multi-scale helix filter */
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
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc) || defined (__sgi)
#include <alloca.h>
#endif
#endif

#include <rsf.h>

#include "createmshelix.h"
#include "createhelix.h"
#include "bound.h"

#include "mshelix.h"
/*^*/

msfilter createmshelix(int ndim    /* number of dimensions */, 
		       int* nd     /* data size [ndim] */, 
		       int* center /* filter center [ndim] */, 
		       int* gap    /* filter gap [ndim] */, 
		       int ns      /* number of scales */, 
		       int *jump   /* filter scaling [ns] */, 
		       int* na     /* filter size [ndim] */)
/*< allocate and output a multiscale helix filter >*/
{
    msfilter msaa;
    sf_filter aa;
    int is, ih, nh, id, n123, nb[SF_MAX_DIM];

    aa = createhelix(ndim, nd, center, gap, na);
    nh = aa->nh;

    msaa = msallocate(nh, ns);
    for (is=0; is < ns; is++) {
	for (ih=0; ih < nh; ih++) /* expanded scale lags */
	    msaa->lag[is][ih] = aa->lag[ih]*jump[is]; 
	
    }
    sf_deallocatehelix(aa);

    n123=1;
    for (id=0; id < ndim; id++) n123 *= nd[id];    
    msaa->mis = sf_boolalloc2(n123,ns);

    aa = msaa->one;
    for (is=0; is < ns; is++) { /* for all scales */
	onescale( is, msaa); /* extract a filter */  
	aa->mis = NULL;
	for (id=0; id < ndim; id++) nb[id] = na[id]*jump[is];
	bound(ndim, nd, nd, nb, aa); /* set up its boundaries */
	for (id=0; id < n123; id++) 
	    msaa->mis[is][id] = aa->mis[id];  /* save them */
	free (aa->mis);
    }
 
    return msaa;
}

/* 	$Id: createmshelix.c 2521 2007-02-02 00:25:42Z sfomel $	 */
