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

#include "runge.h"
#include "alloc.h"

static int dim, nt;
static float dt, **k, *yk;

void sf_runge_init(int dim1 /* dimensionality */, 
		   int n1   /* number of ray tracing steps */, 
		   float d1 /* step in time */)
/*< initialize >*/
{
    dim = dim1;
    nt = n1;
    dt = d1;
    
    yk = sf_floatalloc(dim);
    k = sf_floatalloc2(dim,3);
}

void sf_runge_close(void)
/*< free allocated storage >*/
{
    free(yk);
    free(*k);
    free(k);
}

float sf_ode23 (float t /* time integration */,
		float* tol /* error tolerance */,
		float* y   /* [dim] solution */, 
		void* par  /* parameters for function evaluation */,
		void (*rhs)(void*,float*,float*) 
		/* RHS function */, 
		int (*term)(void*,float*)
	     /* function returning 1 if the ray needs to terminate */)
/*< ODE solver for dy/dt = f where f comes from rhs(par,y,f)
  Note: Value of y is changed inside the function.
  >*/
{
    int i;
    float h, t1;
    bool pass;

    t1 = 0.;
    h = dt;
    
    while (h > 0.) {
	if (t1+h >= t) h=t-t1;

	rhs (par, y, k[0]); 
	for (i=0; i < dim; i++) {
	    yk[i] = y[i] + 0.5*h*k[0][i];
	}
      
	if (term != NULL && term (par, yk)) return t1;
	rhs (par, yk, k[1]); 
	for (i=0; i < dim; i++) {
	    yk[i] = y[i] + h*(2.0*k[1][i]-k[0][i]);
	}
      
	if (term != NULL && term (par, yk)) return t1;
	rhs (par, yk, k[2]); 
	pass = true;
	for (i=0; i < dim; i++) {	    
	    yk[i] = h*(k[0][i]+4.*k[1][i]+k[2][i])/6.0;
	    if (fabsf(yk[i]-h*k[1][i]) > tol[i]) {
		pass = false;
		break;
	    }
	}
	if (pass) {
	    for (i=0; i < dim; i++) {
		y[i] += yk[i];
	    }
	    t1 += h;	
	    if (term != NULL && term (par, y)) return t1;
	} else {
	    h *= 0.5;
	}
    }
  
    return t1;
}


int sf_ode23_step (float* y    /* [dim] solution */, 
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
 
    if (traj != NULL) {
	for (i=0; i < dim; i++) {
	    traj[0][i] = y[i];
	}
    }

    for (it = 0; it < nt; it++) {
	rhs (par, y, k[0]); 
	for (i=0; i < dim; i++) {
	    yk[i] = y[i] + 0.5*dt*k[0][i];
	}
      
	rhs (par, yk, k[1]); 
	for (i=0; i < dim; i++) {
	    yk[i] = y[i] + dt*(2.0*k[1][i]-k[0][i]);
	}
      
	rhs (par, yk, k[2]); 
	for (i=0; i < dim; i++) {
	    y[i] += dt*(k[0][i]+4.*k[1][i]+k[2][i])/6.0;
	    if (traj != NULL) traj[it+1][i] = y[i];
	}
	
	if (term != NULL && term (par, y)) return (it+1);
    }
  
    return 0;
}

/* 	$Id: runge.c 746 2004-08-18 07:38:00Z fomels $	 */
