/* Symplectic Runge-Kutta ray tracing. */
/*
  Copyright (C) 2011 University of Texas at Austin
  
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
/* 
   magic coefficients from
   R.I.McLachlan and P.Atela "The accuracy of symplectic integrators".
   Nonlinearity, 5 (1992), 541-562.
*/

int atela_step (int dim      /* dimensionality */, 
		int nt       /* number of ray tracing steps */, 
		float dt     /* step in time */, 
		float shift  /* complex source shift */,
		bool intime  /* if step in time (not sigma) */, 
		float* x     /* [dim] point location */, 
		float* p     /* [dim] ray parameter vector */,
		void* par     
		/* parameters for slowness functions */,
		void* devz     
		/* parameters for slowness derivative functions */,
		void* devy     
		/* parameters for slowness derivative functions */,
		void (*vgrad)(void*,float*,float*) 
		/* function returning 1/2*(gradient of slowness squared) */, 
		float (*slow2)(void*,float*)
		/* function returning slowness squared */, 
		int (*term)(void*,float*)
		/* function returning 1 if the ray needs to terminate */, 
		float** traj   
		/* [nt+1][dim] - ray trajectory (output) */,
		float** dire   
		/* [nt+1][dim] - ray direction (output) */,
		sf_complex*** dynaM 
		/* [nt+1][dim][dim] - dynamic ray (output) */,
		sf_complex*** dynaN 
		/* [nt+1][dim][dim] - dynamic ray (output) */) 
/*< ray tracing 
  Note:
  1. Values of x and p are changed inside the function.
  2. The trajectory traj is stored as follows:
  {z0,y0,z1,y1,z2,y2,...} in 2-D
  {z0,y0,x0,z1,y1,x1,...} in 3-D (NOTE: disabled)
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
    float f[3], qz[3], qy[3];
    float v, h=0.;
    double sum, one;
    sf_complex M[2][2], N[2][2];

    /* initialize complex matrix */
    v = slow2(par, x);

    M[0][0] = sf_cmplx(0.,sqrtf(v)); M[0][1] = sf_cmplx(0.,0.); M[1][0] = sf_cmplx(0.,0.); M[1][1] = sf_cmplx(0.,sqrtf(v));
    N[0][0] = sf_cmplx(shift,0.);    N[0][1] = sf_cmplx(0.,0.); N[1][0] = sf_cmplx(0.,0.); N[1][1] = sf_cmplx(shift,0.);
 
    /* write initial values */
    if (traj != NULL) {
	for (i=0; i < dim; i++)
	    traj[0][i] = x[i];
    }
    if (dire != NULL) {
	for (i=0; i < dim; i++)
	    dire[0][i] = p[i];
    }
    if (dynaM != NULL) {
	dynaM[0][0][0] = M[0][0]; dynaM[0][0][1] = M[0][1]; dynaM[0][1][0] = M[1][0]; dynaM[0][1][1] = M[1][1];
    }
    if (dynaN != NULL) {
	dynaN[0][0][0] = N[0][0]; dynaN[0][0][1] = N[0][1]; dynaN[0][1][0] = N[1][0]; dynaN[0][1][1] = N[1][1];
    }
    
    /* increment in sigma */
    if (!intime) h = dt;

    for (it = 0; it < nt; it++) {
	/* loop over time steps */

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

	/* symplectic Runge-Kutta */
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
	    if (dire != NULL) dire[it+1][i] = p[i];
	}
     
	/* extrapolate complex matrix */
	/* one step, first-order accuracy */
	vgrad (devz,x,qz);
	vgrad (devy,x,qy);

	/* enforce symmetry grad_x(grad_z(s*s))=grad_z(grad_x(s*s)) */
	qz[1] = qy[0] = 0.5*(qz[1]+qy[0]);

	/* factor 2 in front to compensate for vgrad */
	M[0][0] += 2.*h*(qz[0]*N[0][0]+qz[1]*N[0][1]);
	M[0][1] += 2.*h*(qy[0]*N[0][0]+qy[1]*N[0][1]);
	M[1][0] += 2.*h*(qz[0]*N[1][0]+qz[1]*N[1][0]);
	M[1][1] += 2.*h*(qy[0]*N[1][0]+qy[1]*N[1][1]);

	N[0][0] += h*M[0][0];
	N[0][1] += h*M[0][1];
	N[1][0] += h*M[1][0];
	N[1][1] += h*M[1][1];

	if (dynaM != NULL) {
	    dynaM[it+1][0][0] = M[0][0]; dynaM[it+1][0][1] = M[0][1]; dynaM[it+1][1][0] = M[1][0]; dynaM[it+1][1][1] = M[1][1];
	}
	if (dynaN != NULL) {
	    dynaN[it+1][0][0] = N[0][0]; dynaN[it+1][0][1] = N[0][1]; dynaN[it+1][1][0] = N[1][0]; dynaN[it+1][1][1] = N[1][1];
	}

	/* early termination if reaches domain boundaries */
	if (term != NULL && term (par, x)) return (it+1);
    }
  
    return 0;
}
