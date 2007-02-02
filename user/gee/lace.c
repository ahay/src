/* Estimating prediction-error filter beyond aliasing */
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

#include "lace.h"
#include "createhelix.h" 
#include "bound.h"
#include "pef.h"

sf_filter lace_pef(int dim     /* number of dimensions */, 
		   float *dd   /* data */, 
		   int jump    /* filter stretch */, 
		   int n       /* data size */, 
		   int *nd     /* data dimensions [dim] */, 
		   int *center /* filter center [dim] */, 
		   int *gap    /* filter gap [dim] */, 
		   int *na     /* filter size [dim] */)  
/*< estimate PEF >*/
{
    int *savelags, ii[SF_MAX_DIM]; /* holding place */
    int ih, nh, lag0, j;
    sf_filter aa;

    aa = createhelix(dim, nd, center, gap, na);  
    savelags = aa->lag;
    nh = aa->nh;

    aa->lag = sf_intalloc(nh); /* prepare interlaced helix */
    lag0 = sf_cart2line(dim, na, center);

    for (ih=0; ih < nh; ih++) {	/* sweep through the filter */
	sf_line2cart(dim, na, ih+lag0+1, ii);
	for (j=0; j < dim; j++) {
	    ii[j] -= center[j];
	}
	ii[0] *= jump;  /* interlace on 1-axis */
	aa->lag[ih] = sf_cart2line(dim, nd, ii);
    }
    na[0] *= jump;
    bound(dim, nd, nd, na, aa);  /* define  aa->mis */
    na[0] /= jump;
    
    find_pef(n, dd, aa, nh*2);    /* estimate aa coefficients */
    free(aa->lag);   
    aa->lag = savelags;		  /* restore filter lags */

    return aa;
}

