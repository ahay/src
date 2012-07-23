/* Ray tracing interface. */
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

#include "raytrace.h"
#include "grid2.h"
#include "grid3.h"
#include "atela.h"

#ifndef _raytrace_h

typedef struct RayTrace* raytrace;
/* abstract data type */
/*^*/

#endif

struct RayTrace {
    bool sym;
    int dim, nt;
    float dt, z0;
    grid3 grd3;
    grid2 grd2;
};
/* concrete data type */

static void iso_rhs(void* par, float* y, float* f)
/* right-hand side for isotropic raytracing */
{    
    raytrace rt;
    int i, dim;
    float s2, sds[3];
	
    rt = (raytrace) par;
    dim = rt->dim;
	
    switch (dim) {
		case 2:
			s2 = grid2_vel(rt->grd2,y);
			grid2_vgrad(rt->grd2,y,sds);
			break;
		case 3:
			s2 = grid3_vel(rt->grd3,y);
			grid3_vgrad(rt->grd3,y,sds);
			break;
		default:
			s2 = 0.;
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }
	
    for (i=0; i < dim; i++) {
		f[i]   = y[i+dim]/s2; /* p/s^2 */
		f[i+dim] = sds[i]/s2; /* 1/2 grad(s^2)/s^2 */
    }
}

static int term(void* par, float* y)
/* grid termination */
{
    raytrace rt;
    
    rt = (raytrace) par;
	
    switch (rt->dim) {
		case 2:
			return grid2_term(rt->grd2,y);
		case 3:
			return grid3_term(rt->grd3,y);
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,rt->dim);
			return 0;
    }
}

raytrace raytrace_init(int dim            /* dimensionality (2 or 3) */, 
					   bool sym,          /* if symplectic */
					   int nt             /* number of ray tracing steps */, 
					   float dt           /* ray tracing step (in time) */,
					   int* n             /* slowness dimensions [dim] */, 
					   float* o, float* d /* slowness grid [dim] */,
					   float* slow2       /* slowness squared [n3*n2*n1] */, 
					   int order          /* interpolation order */)
/*< Initialize ray tracing object. 
 * Increasing order increases accuracy but
 decreases efficiency. Recommended values: 3 or 4.
 * slow2 can be changed or deallocated after
 raytrace_init.
 >*/
{
    raytrace rt;
    
    rt = (raytrace) sf_alloc (1,sizeof(*rt));
    
    rt->dim = dim;
    rt->sym = sym;
    rt->nt = nt;
    rt->dt = dt;
    rt->z0 = o[0];
    
    switch (dim) {
		case 2:
			rt->grd2 = grid2_init (n[0], o[0], d[0], 
								   n[1], o[1], d[1],
								   slow2, order);
			break;
		case 3:
			rt->grd3 = grid3_init (n[0], o[0], d[0], 
								   n[1], o[1], d[1],
								   n[2], o[2], d[2],
								   slow2, order);
			break;
		default:
			sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }
	
    return rt;
}

void raytrace_close (raytrace rt)
/*< Free internal storage >*/
{
    switch (rt->dim) {
		case 2:
			grid2_close (rt->grd2);
			break;
		case 3:
			grid3_close (rt->grd3);
			break;
    }
    free (rt);
}

int trace_ray (raytrace rt  /* ray tracing object */, 
			   float* x     /* point location {z,y,x} [dim] */, 
			   float* p     /* ray parameter vector [dim] */, 
			   float** traj /* output ray trajectory [nt+1,dim] */)
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
    float y[6], s2;
	
    dim = rt->dim;
    nt = rt->nt;
	
    if (!rt->sym) {
		switch (dim) {
			case 2:
				s2 = grid2_vel(rt->grd2,x);
				break;
			case 3:
				s2 = grid3_vel(rt->grd3,x);
				break;
			default:
				s2 = 0.;
				sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
		}
		
		for (i=0; i < dim; i++) {
			y[i] = x[i];
			y[i+dim] = p[i]*sqrtf(s2);
		}
		
		sf_runge_init(2*dim, nt, rt->dt);
		it = sf_ode23_step (y, rt,iso_rhs,term,traj);
		sf_runge_close();
		
		for (i=0; i < dim; i++) {
			x[i] = y[i];
			p[i] = y[i+dim];
		}
    } else {
		switch (dim) {
			case 2:
				it = atela_step (dim, nt, rt->dt, true, x, p, 
								 rt->grd2, 
								 grid2_vgrad, grid2_vel, grid2_term, traj);
				break;
			case 3:
				it = atela_step (rt->dim, nt, rt->dt, true, x, p, 
								 rt->grd3, 
								 grid3_vgrad, grid3_vel, grid3_term, traj);
				break;
			default:
				sf_error("%s: cannot handle %d dimensions",__FILE__,rt->dim);
				break;
		}
    }
    
    if (it > 0 && x[0] > rt->z0) {
		return (-it); /* exit through the side or bottom */
    } else {
		return it;
    }
}

/* 	$Id$	 */
