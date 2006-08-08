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

void sf_tinysolver (sf_operator oper   /* linear operator */, 
		    sf_solverstep solv /* stepping function */, 
		    int nx             /* size of x */, 
		    int ny             /* size of dat */, 
		    float* x           /* estimated model */, 
		    const float* dat   /* data */, 
		    int niter          /* number of iterations */)
/*< Generic linear solver. Solves oper{x} =~ dat >*/
{
    int i, iter;
    float *g, *rr, *gg;
    bool forget = false;

    g =  sf_floatalloc (nx);
    rr = sf_floatalloc (ny);
    gg = sf_floatalloc (ny);

    for (i=0; i < ny; i++) rr[i] = - dat[i];
    for (i=0; i < nx; i++) x[i] = 0.0;

    for (iter=0; iter < niter; iter++) {
	oper (true, false, nx, ny, g, rr);
	oper (false, false, nx, ny, g, gg);

	solv (forget, nx, ny, x, g, rr, gg);
    }

    free (g);
    free (rr);
    free (gg);
}
