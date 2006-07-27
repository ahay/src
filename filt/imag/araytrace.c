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
#include "runge.h"

#ifndef _araytrace_h

typedef struct aRayTrace* araytrace;
/* abstract data type */
/*^*/

#endif

struct aRayTrace {
    int dim, nt;
    float dt, z0;
    grid2a grd2;
};
/* concrete data type */

static void aniso_rhs(void* par, float* y, float* f)
/* right-hand side for anisotropic raytracing */
{    
    araytrace rt;
    int dim;

    rt = (araytrace) par;
    dim = rt->dim;

    switch (dim) {
	case 2:
	    grid2a_rhs(rt->grd2,y,f);
	    break;
	default:
	    sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }
}

static int term(void* par, float* y)
/* grid termination */
{
    araytrace rt;
    
    rt = (araytrace) par;

    switch (rt->dim) {
	case 2:
	    return grid2a_term(rt->grd2,y);
	default:
	    sf_error("%s: Cannot raytrace with dim=%d",__FILE__,rt->dim);
	    return 0;
    }
}

araytrace araytrace_init(int dim            /* dimensionality (2 or 3) */, 
			 int nt             /* number of ray tracing steps */, 
			 float dt           /* ray tracing step (in time) */,
			 int* n             /* velocity dimensions [dim] */, 
			 float* o, float* d /* velocity grid [dim] */,
			 float* vz2         /* vertical velocity squared [n3*n2*n1] */, 
			 float* vx2         /* horizontal velocity squared [n3*n2*n1] */,
			 float* q           /* analepticity */,
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
    
    switch (dim) {
	case 2:
	    rt->grd2 = grid2a_init (n[0], o[0], d[0], 
				    n[1], o[1], d[1],
				    vz2, vx2, q, order);
	    break;
	default:
	    sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }

    return rt;
}

void araytrace_close (araytrace rt)
/*< Free internal storage >*/
{
    switch (rt->dim) {
	case 2:
	    grid2a_close (rt->grd2);
	    break;
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

    runge_init(2*dim, nt, rt->dt);
    it = ode23_step (y, rt, aniso_rhs, term, traj);
    runge_close();
    
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

/* 	$Id: raytracing.c 775 2004-09-03 19:44:45Z fomels $	 */
