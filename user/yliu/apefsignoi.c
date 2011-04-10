/* Signal and noise separation with adaptive prediction-error filters */
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

#include "apefsignoi.h"
#include "apefs.h" 

void apefsignoi_lop(int niter         /* number of iterations */, 
		    int ssf1, int ssf2,
		    int nnf1, int nnf2,
		    int nn1, int nn2,
		    float *ssfilt,
		    float *nnfilt,
		    float *xx         /* model */, 
		    float *yy         /* data  */,
		    float eps         /* regularization parameter */,
		    bool verb         /* verbosity flag */) 
/*< linear operator >*/
{
    int ix;
    float *dd;

    /*  regularized */
    dd = sf_floatalloc(nn1*nn2);
    for (ix=0; ix < nn1*nn2; ix++) {
	dd[ix] = 0.;
    }
    
    apefs_init(ssfilt, nnfilt, ssf1, ssf2, nn1, nn2, nnf1, nnf2);

    npefs_lop(false,false,nn1*nn2,nn1*nn2,yy,dd);
    /* N^2{d} */

    sf_solver_reg (npefs_lop, sf_cgstep, spefs_lop, nn1*nn2, nn1*nn2, nn1*nn2, 
		   xx, dd, niter, eps, "verb", verb, "end");
    /* N^2{x} =~ N^2{d} && eps S{x} =~ 0 */
    free(dd);
    sf_cgstep_close();
}

/* 	$Id$	 */

