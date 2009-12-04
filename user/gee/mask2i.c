/* Missing data interpolation with one or two PEFs */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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

#include "heliarr.h"

void maski (int niter       /* number of iterations */, 
	    int nx          /* data size */, 
	    const float *dd /* data */, 
	    float *xx       /* model */, 
	    int *known      /* mask for known data */, 
	    sf_filter aa1   /* first PEF */, 
	    sf_filter aa2   /* second PEF */) 
/*< interpolate >*/
{
    bool *mask;
    int ix;

    mask = sf_boolalloc(nx);
    for (ix=0; ix < nx; ix++) {
	mask[ix] = (bool) known[ix];
    }

    sf_mask_init(mask);

    if (NULL != aa2) {
	heliarr_init (aa1, aa2);
	sf_solver_reg(sf_mask_lop,sf_cgstep,
		      heliarr_lop,2*nx,nx,nx,xx,dd,niter,1.,"end");
    } else {
	sf_helicon_init (aa1);
	sf_solver_reg(sf_mask_lop,sf_cgstep,
		      sf_helicon_lop,nx,nx,nx,xx,dd,niter,1.,"end");
    }

    sf_cgstep_close();
    free(mask);
}





