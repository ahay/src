/* Prediction-error filter estimation */
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

#include "pef.h"
#include "hconest.h"

void find_pef(int nd       /* data size */, 
	      float* dd    /* data [nd] */, 
	      sf_filter aa /* estimated filter */, 
	      int niter    /* number of iterations */)
/*< find PEF >*/ 
{
    hconest_init( dd, aa);
    sf_solver(hconest_lop, sf_cgstep, aa->nh, nd, aa->flt, dd, 
	      niter, "x0", aa->flt, "end");
    sf_cgstep_close();
}

