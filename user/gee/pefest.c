/* Estimate a PEF avoiding zeros and bursty noise on input. */
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

#include "misinput.h"
#include "pef.h"

void pefest(int niter    /* number of iterations */,
	    int ny       /* data size */,
	    float *yy    /* input data */, 
	    sf_filter aa /* estimated PEF */) 
/*< estimate PEF avoiding zeros and bursty noise >*/
{
    int iy, *mask;
    float *rr, *rabs, rbar;

    rr = sf_floatalloc(ny);
    rabs = sf_floatalloc(ny);
    mask = sf_intalloc(ny);

    sf_helicon_init(aa);                /* starting guess */
    sf_helicon_lop(false,false,ny,ny,yy,rr);
    for (iy=0; iy < ny; iy++) {
	rabs[iy] = fabsf(rr[iy]);
    }
    rbar = sf_quantile(ny/3,ny,rabs); /* rbar = (r safe below rbar) */
    for (iy=0; iy < ny; iy++) {
	if (aa->mis[iy]) {
	    yy[iy] = 0.;
	    mask[iy] = 0;
	} else {
	    mask[iy] = fabsf(rr[iy]) < 5.*rbar;
	}
    }	
    find_mask(ny,mask,aa);
    find_pef(ny,yy,aa,niter); 
 
    free(rr);
    free(rabs);
    free(mask);
}
