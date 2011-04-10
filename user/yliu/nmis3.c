/* 3-D missing data interpolation with nonstationary filters */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "nmis3.h"
#include "mmmult3.h" 

void nmis3(int niter         /* number of iterations */, 
	   int nf1, int nf2, 
	   int nf3, int nf4,
	   int nf5, int nf6,
	   float *filt,
	   float *xx         /* model */, 
	   const bool *known /* mask for known data */,
	   float eps         /* regularization parameter */,
	   bool verb         /* verbosity flag */) 
/*< 3-D interpolation >*/
{
    int ix;
    float *dd;

    /*  regularized */
    dd = sf_floatalloc(nf4*nf5*nf6);
    for (ix=0; ix < nf4*nf5*nf6; ix++) {
	dd[ix]=0.;
    }
    
    mmmult3_init(filt, nf1, nf2, nf3, nf4, nf5, nf6);
    sf_solver (mmmult3_lop, sf_cgstep, nf4*nf5*nf6, nf4*nf5*nf6, 
	       xx, dd, niter, "known", known, "x0", xx, "verb", verb, "end");
    free(dd);
    sf_cgstep_close();
    
}

/* 	$Id$	 */

