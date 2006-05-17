/* Runge-Kutta ODE solvers. */
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

#include <stdlib.h>
#include <math.h>

#include <rsf.h>
/*^*/

#include "runge.h"

int ode23_step (int dim     /* dimensionality */, 
		int nt      /* number of ray tracing steps */, 
		float dt    /* step in time */, 
		float* y    /* [dim] solution */, 
		void* par   /* parameters for function evaluation */,
		void (*rhs)(void*,float*,float*) 
		/* RHS function */, 
		int (*term)(void*,float*)
		/* function returning 1 if the ray needs to terminate */, 
		float** traj /* [nt+1][dim] - ray trajectory (output) */) 
/*< ODE solver for dy/dt = f where f comes from rhs(par,y,f)
  Note:
  1. Value of y is changed inside the function.
  2. The output code for it = ode23_step(...)
  it=0 - ray traced to the end without termination
  it>0 - ray terminated
  The total traveltime along the ray is 
  nt*dt if (it = 0); it*dt otherwise 
  >*/
{
    int it, i;
    float f[3], g[3], gi;
 
    if (traj != NULL) {
	for (i=0; i < dim; i++) {
	    traj[0][i] = y[i];
	}
    }

    for (it = 0; it < nt; it++) {
	rhs (par, y, f); /* k1 */
	for (i=0; i < dim; i++) {
	    gi = 0.5*dt*f[i];
	    y[i] += gi;
	    g[i] = gi;
	}
      
	rhs (par, y, f); /* k2 */
	for (i=0; i < dim; i++) {
	    gi = 0.75*dt*f[i];
	    y[i] += gi;
	    g[i] += gi;
	}
      
	rhs (par, y, f); /* k3 */
	for (i=0; i < dim; i++) {
	    y[i] += g[i] + 4.*(g[i]+dt*f[i])/9.;
	    if (traj != NULL) traj[it+1][i] = y[i];
	}
     
	/* Error: -5/72 k1 + 1/12 k2 + 1/9 k3 - 1.8 k4 */

	if (term != NULL && term (par, y)) return (it+1);
    }
  
    return 0;
}

/* 	$Id: atela.c 746 2004-08-18 07:38:00Z fomels $	 */
