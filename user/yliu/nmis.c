/* Missing data interpolation with nonstationary filters */
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

#include "nmis.h"
#include "mmmult.h" 

void nmis(int niter         /* number of iterations */, 
	  int nf1, int nf2, 
	  int nf3, int nf4,
	  float *filt,
	  float *xx         /* model */, 
	  const bool *known /* mask for known data */,
	  float eps         /* regularization parameter */,
	  bool verb         /* verbosity flag */) 
/*< interpolate >*/
{
    int ix;
    float *dd;

    /*  regularized */
    dd = sf_floatalloc(nf3*nf4);
    for (ix=0; ix < nf3*nf4; ix++) {
	dd[ix]=0.;
    }
    
    mmmult_init(filt, nf1, nf2, nf3, nf4);
    sf_solver (mmmult_lop, sf_cgstep, nf3*nf4, nf3*nf4, xx, dd, niter, 
	       "known", known, "x0", xx, "verb", verb, "end");
    free(dd);
    sf_cgstep_close();

}


