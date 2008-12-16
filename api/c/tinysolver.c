/* Solver function for iterative least-squares optimization. */
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
#include "_solver.h"
#include "alloc.h"
#include "tinysolver.h"

void sf_tinysolver (sf_operator Fop       /* linear operator */, 
		    sf_solverstep stepper /* stepping function */, 
		    int nm                /* size of model */, 
		    int nd                /* size of data */, 
		    float* m              /* estimated model */,
		    const float* m0       /* starting model */,
		    const float* d        /* data */, 
		    int niter             /* iterations */)
/*< Generic linear solver. Solves oper{x} =~ dat >*/
{
    int i, iter;
    float *g, *rr, *gg;

    g =  sf_floatalloc (nm);
    rr = sf_floatalloc (nd);
    gg = sf_floatalloc (nd);

    for (i=0; i < nd; i++) rr[i] = - d[i];
    if (NULL==m0) {
	for (i=0; i < nm; i++) m[i] = 0.0;
    } else {
	for (i=0; i < nm; i++) m[i] = m0[i];
	Fop (false, true, nm, nd, m, rr);
    }
    
    for (iter=0; iter < niter; iter++) {
	Fop (true, false, nm, nd, g, rr);
	Fop (false, false, nm, nd, g, gg);
	
	stepper (false, nm, nd, m, g, rr, gg);
    }
    
    free (g);
    free (rr);
    free (gg);
}
