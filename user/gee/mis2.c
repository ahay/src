/* Missing data interpolation with helical filters */
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

#include "mis2.h"
#include "helicon.h" 
#include "polydiv.h" 

void mis2(int niter         /* number of iterations */, 
	  int nx            /* model size */, 
	  float *xx         /* model */, 
	  sf_filter aa      /* helix filter */, 
	  const bool *known /* mask for known data */,
	  float eps         /* regularization parameter */,
	  bool doprec       /* to apply preconditioning */) 
/*< interpolate >*/
{
    int ix;
    float *dd;

    if (doprec) {                          /*  preconditioned */
	sf_mask_init(known);
	polydiv_init(nx, aa);
	sf_solver_prec(sf_mask_lop, sf_cgstep, polydiv_lop, 
		       nx, nx, nx, xx, xx, niter, eps, "end");
	polydiv_close();
    } else {                               /*  regularized */
	dd = sf_floatalloc(nx);
	for (ix=0; ix < nx; ix++) {
	    dd[ix]=0.;
	}

	helicon_init(aa);
	sf_solver (helicon_lop, sf_cgstep, nx, nx, xx, dd, niter, 
		   "known", known, "x0", xx, "end");
	free(dd);
    }
    sf_cgstep_close();
}

/* 	$Id$	 */

