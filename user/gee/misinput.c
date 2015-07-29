/* Create a missing data mask for a helical filter */
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

#include "misinput.h" 

void find_mask(int n            /* data size */, 
	       const int *known /* mask for known data [n] */, 
	       sf_filter aa     /* helical filter */) 
/*< create a filter mask >*/
{
    int i, ih;
    float *rr, *dfre;

    rr = sf_floatalloc(n);
    dfre = sf_floatalloc(n);

    for (i=0; i < n; i++) {
	dfre[i] = known[i]? 0.:1.;
    }
    
    sf_helicon_init(aa);

    for (ih=0; ih < aa->nh; ih++) {
	aa->flt[ih] = 1.;
    }

    sf_helicon_lop(false,false,n,n,dfre,rr);

    for (ih=0; ih < aa->nh; ih++) {
	aa->flt[ih] = 0.;
    }

    for (i=0; i < n; i++) {
	if ( rr[i] > 0.) aa->mis[i] = true;	
    }

    free (rr);
    free (dfre);
}

/* 	$Id: misinput.c 7107 2011-04-10 02:04:14Z ivlad $	 */
