/* INVerse INTerpolation in 1-D. */
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

#include "lint1.h"
#include "tcai1.h"

void invint1(int niter                  /* number of iterations */,
	     int nd                     /* data size */,
	     float *coord               /* data coordinates */, 
	     const float *dd            /* data values */, 
	     int n1, float o1, float d1 /* model grid */, 
	     int na, const float *aa    /* filter */, 
	     float *mm                  /* estimated model */, 
	     float eps                  /* regularization */)
/*< inverse interpolation >*/
{
    lint1_init( o1, d1, coord); /* interpolation */
    tcai1_init( na, aa);         /* filtering */

    sf_solver_reg(lint1_lop, sf_cgstep, tcai1_lop, 
		  n1+na, n1, nd, mm, dd, niter, eps);
    sf_cgstep_close( );
}
