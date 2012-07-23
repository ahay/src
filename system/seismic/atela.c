/* Symplectic Runge-Kutta ray tracing. */
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

#include "atela.h"

static const double a1 = 0.5153528374311229364;
static const double a2 = -0.085782019412973646;
static const double a3 = 0.4415830236164665242;
static const double a4 = 0.1288461583653841854;
static const double b1 = 0.1344961992774310892;
static const double b2 = -0.2248198030794208058;
static const double b3 = 0.7563200005156682911;
static const double b4 = 0.3340036032863214255;
/* magic coefficients from
 R.I.McLachlan and P.Atela "The accuracy of symplectic integrators".
 Nonlinearity, 5 (1992), 541-562.
 */

int atela_step (int dim     /* dimensionality */, 
				int nt      /* number of ray tracing steps */, 
				float dt    /* step in time */, 
				bool intime /* if step in time (not sigma) */, 
				float* x    /* [dim] point location */, 
				float* p    /* [dim] ray parameter vector */,
				void* par   /* parameters for slowness functions */,
				void (*vgrad)(void*,float*,float*) 
				/* function returning 1/2*(gradient of slowness squared) */, 
				float (*slow2)(void*,float*)
				/* function returning slowness squared */, 
				int (*term)(void*,float*)
				/* function returning 1 if the ray needs to terminate */, 
				float** traj /* [nt+1][dim] - ray trajectory (output) */) 
/*< ray tracing 
 Note:
 1. Values of x and p are changed inside the function.
 2. The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 3. Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 4. The output code for it = atela_step(...)
 it=0 - ray traced to the end without termination
 it>0 - ray terminated
 The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int it, i;
    float f[3];
    float v, h=0.;
    double sum, one;
	
    if (traj != NULL) {
		for (i=0; i < dim; i++) {
			traj[0][i] = x[i];
		}
    }
	
    if (!intime) h = dt;
	
    for (it = 0; it < nt; it++) {
		v = slow2(par, x);
		sum = 0.;
		for (i=0; i < dim; i++) {
			sum += p[i]*p[i];
		}
		one = sqrt(v/sum);
		for (i=0; i < dim; i++) {
			p[i] *= one; /* enforce the eikonal equation */
		}
		
		if (intime) h = dt/v;
		
		vgrad (par, x, f);
		for (i=0; i < dim; i++) {
			p[i] += b1*h*f[i];
			x[i] += a1*h*p[i];
		}
		
		vgrad (par, x, f);
		for (i=0; i < dim; i++) {
			p[i] += b2*h*f[i];
			x[i] += a2*h*p[i];
		}
		
		vgrad (par, x, f);
		for (i=0; i < dim; i++) {
			p[i] += b3*h*f[i];
			x[i] += a3*h*p[i];
		}
		
		vgrad (par, x, f);
		for (i=0; i < dim; i++) {
			p[i] += b4*h*f[i];
			x[i] += a4*h*p[i];
			if (traj != NULL) traj[it+1][i] = x[i];
		}
		
		if (term != NULL && term (par, x)) return (it+1);
    }
	
    return 0;
}

/* 	$Id$	 */
