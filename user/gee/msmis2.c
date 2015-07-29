/* Multi-scale missing data interpolation */
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
/*^*/

#include "msmis2.h"
#include "mshelicon.h" 

#include "mshelix.h"
/*^*/

void msmis2(int niter         /* number of iterations */, 
	    int nx            /* data size */, 
	    int ns            /* number of scales */, 
	    float *xx         /* data */, 
	    msfilter aa       /* filter */, 
	    const bool *known /* mask for known data */)
/*< interpolate >*/ 
{
    int ix, nxs;
    float *dd;

    nxs = ns*nx;

    dd = sf_floatalloc(nxs);
    for (ix=0; ix < nxs; ix++) {
	dd[ix]=0.;
    }
    
    mshelicon_init(aa);
    sf_solver (mshelicon_lop, sf_cgstep, nx, nxs, xx, dd, niter, 
	       "known", known, "x0", xx, "end");
    
    sf_cgstep_close();
    free(dd);
}

/* 	$Id: msmis2.c 841 2004-10-25 13:08:40Z fomels $	 */
