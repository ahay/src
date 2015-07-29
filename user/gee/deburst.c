/* Remove burst noise, 1-D */
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

#include "deburst.h"
#include "icai1.h"

void deburst (int n           /* data length */, 
	      int niter       /* number of iterations */, 
	      sf_weight wght  /* weight operator */, 
	      float eps       /* regularization */, 
	      const float *dd /* input data */, 
	      float *hh       /* output model */) 
/*< debursting by optimization >*/
{
    float aa[] = {-1.,2.,-1.}; /* laplacian filter */

    icai1_init (3, aa, 1); 
    sf_solver_reg(sf_copy_lop, sf_cgstep, icai1_lop, n, n, n, hh, dd, 
		  niter, eps, "wght", wght, "nfreq", 1, "nmem", 0, "end");
    sf_cgstep_close();
}
