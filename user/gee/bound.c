/* Helix filter boundary conditions. */
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
#include <stdio.h>

#include <rsf.h>

#include "bound.h"
#include "regrid.h"

void bound (int dim         /* number of dimensions */, 
	    const int *nold /* old data coordinates [dim] */, 
	    const int *nd   /* new data coordinates [dim] */, 
	    const int *na   /* filter box size [dim] */, 
	    const sf_filter aa /* helix filter */) 
/*< Mark helix filter outputs where input is off data. >*/
{
    int iy, my, ib, mb, i, nb[SF_MAX_DIM], ii[SF_MAX_DIM];
    float *xx, *yy;
    
    my = mb = 1;
    for (i=0; i < dim; i++) {
	nb[i] = nd[i] + 2*na[i]; /* nb is a bigger space. */   
	mb *= nb[i];
	my *= nd[i];
    }

    xx = sf_floatalloc(mb); yy = sf_floatalloc(mb);

    for (ib=0; ib < mb; ib++) { 
	sf_line2cart(dim, nb, ib, ii);    
	xx[ib] = 0.;  	
	for (i=0; i < dim; i++) 
	    if(ii[i]+1 <= na[i] ||  ii[i]+1 > nb[i]-na[i]) {
		xx[ib] = 1.; 
		break;
	    }
    }
    sf_helicon_init( aa);		  
    regrid(dim, nold, nb, aa);  
    for (i=0; i < aa->nh; i++) aa->flt[i] = 1.;		
    /* apply filter */
    sf_helicon_lop(false, false, mb, mb, xx, yy); 
    regrid(dim, nb, nd, aa);  
    for (i=0; i < aa->nh; i++) aa->flt[i] = 0.;

    aa->mis = sf_boolalloc(my); /* attach missing designation */
    for (iy = 0; iy < my; iy++) {  /* map to padded space */
	sf_line2cart(dim, nd, iy, ii);
	for (i=0; i < dim; i++) ii[i] += na[i];
	ib = sf_cart2line(dim, nb, ii);
	aa->mis[iy] = (bool) (yy[ib] > 0.);  
    }
    
    free (xx); free (yy);
} 

/* 	$Id: bound.c 3860 2008-11-06 15:23:21Z sfomel $	 */
