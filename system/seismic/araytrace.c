/* Ray tracing interface for VTI. */
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

#include "araytrace.h"
#include "grid2a.h"
#include "grid3a.h"

#ifndef _araytrace_h

typedef struct aRayTrace* araytrace;
/* abstract data type */
/*^*/

#endif

struct aRayTrace {
    int dim, nt;
    float dt, z0;
    grid2a grd2;
    grid3a grd3;
};
/* concrete data type */

static void aniso_rhs(void* par, float* y, float* f)
/* right-hand side for anisotropic raytracing */
{    
    araytrace rt;
    int dim;

    rt = (araytrace) par;
    dim = rt->dim;

    if (2==dim) {
	grid2a_rhs(rt->grd2,y,f);
    } else {
	grid3a_rhs(rt->grd3,y,f);
    }
}

static int term(void* par, float* y)
/* grid termination */
{
    araytrace rt;
    int dim;

    rt = (araytrace) par;
    dim = rt->dim;

    if (2==dim) {
	return grid2a_term(rt->grd2,y);
    } else {
	return grid3a_term(rt->grd3,y);
    }
}

araytrace araytrace_init(int dim            /* dimensionality (2 or 3) */, 
			 int nt             /* number of ray tracing steps */, 
			 float dt           /* ray tracing step (in time) */,
			 int* n             /* velocity dimensions [dim] */, 
			 float* o, float* d /* velocity grid [dim] */,
			 float* vz2         /* vertical velocity squared [n3*n2*n1] */, 
			 float* vx2         /* horizontal velocity squared [n3*n2*n1] */,
			 float* q           /* anallepticity [n3*n2*n1] */,
			 int order          /* interpolation order */)
/*< Initialize ray tracing object. 
 * Increasing order increases accuracy but
 decreases efficiency. Recommended values: 3 or 4.
 * vz2, vx2, and q can be changed or deallocated after
 araytrace_init.
 >*/
{
    araytrace rt;
    
    rt = (araytrace) sf_alloc (1,sizeof(*rt));
    
    rt->dim = dim;
    rt->nt = nt;
    rt->dt = dt;
    rt->z0 = o[0];
    
    if (2==dim) {
	rt->grd2 = grid2a_init (n[0], o[0], d[0], 
				n[1], o[1], d[1],
				vz2, vx2, q, order);
	rt->grd3=NULL;
    } else {
	rt->grd2=NULL;
	rt->grd3 = grid3a_init (n[0], o[0], d[0], 
				n[1], o[1], d[1],
				n[2], o[2], d[2],
				vz2, vx2, q, order);
    }

    return rt;
}

void araytrace_close (araytrace rt)
/*< Free internal storage >*/
{
    if (2==rt->dim) {
	grid2a_close (rt->grd2);
    } else {
	grid3a_close (rt->grd3);
    }
    free (rt);
}

int trace_aray (araytrace rt  /* ray tracing object */, 
		float* x      /* point location {z,y,x} [dim] */, 
		float* p      /* ray parameter vector [dim] */, 
		float** traj  /* output ray trajectory [nt+1,dim] */)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6];

    dim = rt->dim;
    nt = rt->nt;

    for (i=0; i < dim; i++) {
	y[i] = x[i];
	y[i+dim] = p[i];
    }

    sf_runge_init(2*dim, nt, rt->dt);
    it = sf_ode23_step (y, rt, aniso_rhs, term, traj);
    sf_runge_close();
    
    for (i=0; i < dim; i++) {
	x[i] = y[i];
	p[i] = y[i+dim];
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}

/* 	$Id: araytrace.c 5023 2009-11-23 01:19:26Z sfomel $	 */
