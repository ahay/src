/* Dynamic ray tracing interface. */
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

#include <rsf.h>

#include "draytrace.h"
#include "grid2.h"
#include "atela.h"

#ifndef _draytrace_h

typedef struct RayTrace* raytrace;
/* abstract data type */
/*^*/

#endif

struct RayTrace {
    /* header information */
    int dim, nt;
    float dt, z0;

    /* ENO grid for s*s */
    grid2 grd2;
    /* ENO grid for s*grad_z(s) */
    grid2 derz;
    /* ENO grid for s*grad_y(s) */
    grid2 dery;
};
/* concrete data type */

raytrace trace_init(int dim            /* dimensionality */, 
		    int nt             /* number of ray tracing steps */, 
		    float dt           /* ray tracing step (in time) */,
		    int* n             /* slowness dimensions [dim] */, 
		    float* o, float* d /* slowness grid [dim] */,
		    float* slow2       /* slowness squared [n0*n1] */, 
		    int order          /* interpolation order */)
/*< Initialize ray tracing object. 
 * Increasing order increases accuracy but
 decreases efficiency. Recommended values: 3 or 4.
 * slow2 can be changed or deallocated after raytrace_init. 
 >*/
{
    raytrace rt;
    float **deriv;
    int iz, iy;
    float temp[2], *sortz, *sorty;

    rt = (raytrace) sf_alloc (1,sizeof(*rt));

    rt->dim = dim;
    rt->nt = nt;
    rt->dt = dt;
    rt->z0 = o[0];
    
    switch (dim) {
	case 2:
	    /* supply s*s to ENO */
	    rt->grd2 = grid2_init (n[0], o[0], d[0], 
				   n[1], o[1], d[1],
				   slow2, order);

	    /* pre-compute s*grad_z(s) and s*grad_y(s) on grid points */
	    deriv = sf_floatalloc2(2,n[0]*n[1]);
	    for (iy=0; iy < n[1]; iy++) {
		temp[1] = o[1]+iy*d[1];

		for (iz=0; iz < n[0]; iz++) {
		    temp[0] = o[0]+iz*d[0];
		    
		    grid2_vgrad ((void *)rt->grd2,temp,deriv[iy*n[0]+iz]);
		}
	    }

	    sortz = sf_floatalloc(n[0]*n[1]);
	    sorty = sf_floatalloc(n[0]*n[1]);
	    for (iz=0; iz < n[0]*n[1]; iz++) {
		sortz[iz] = deriv[iz][0];
		sorty[iz] = deriv[iz][1];
	    }

	    /* supply s*grad_z(s) and s*grad_y(s) to ENO */
	    rt->derz = grid2_init (n[0], o[0], d[0], 
				   n[1], o[1], d[1],
				   sortz, order);
	    rt->dery = grid2_init (n[0], o[0], d[0], 
				   n[1], o[1], d[1],
				   sorty, order);
	    break;
	default:
	    sf_error("%s: Cannot raytrace with dim=%d",__FILE__,dim);
    }

    return rt;
}

void trace_close (raytrace rt)
/*< Free internal storage >*/
{
    switch (rt->dim) {
	case 2:
	    grid2_close (rt->grd2);
	    grid2_close (rt->derz);
	    grid2_close (rt->dery);
	    break;
    }
    free (rt);
}

int trace_ray (raytrace rt  /* ray tracing object */, 
	       float* x     /* point location {z,y,x} [dim] */, 
	       float* p     /* ray parameter vector [dim] */, 
	       float shift  /* complex source shift */,
	       float** traj /* output ray trajectory [nt+1,dim] */,
	       float** dire /* output ray direction [nt+1,dim] */,
	       sf_complex*** dynaM /* output dynamic ray [nt+1,dim,dim] */,
	       sf_complex*** dynaN /* output dynamic ray [nt+1,dim,dim] */)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D (NOTE: disabled)
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.

 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D (NOTE: disabled)
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians

 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int dim, it=0, nt;

    dim = rt->dim;
    nt = rt->nt;
    
    switch (dim) {
	case 2:
	    it = atela_step (dim, nt, rt->dt, shift, true, x, p,
			     rt->grd2, rt->derz, rt->dery,
			     grid2_vgrad, grid2_vel, grid2_term, 
			     traj, dire, dynaM, dynaN);
	    break;
	default:
	    sf_error("%s: cannot handle %d dimensions",__FILE__,rt->dim);
	    break;
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}

void dray_assemble (sf_complex** dynaM, sf_complex** dynaN, sf_complex** dynaK)
/*< assemble dynamic matrix >*/
{
    sf_complex det, inv[2][2];

    /* determinant */
    det = dynaN[0][0]*dynaN[1][1]-dynaN[1][0]*dynaN[0][1];

    /* inv(dynaN) */
    inv[0][0] = dynaN[1][1]/det;
    inv[1][0] = -dynaN[1][0]/det;
    inv[0][1] = -dynaN[0][1]/det;
    inv[1][1] = dynaN[0][0]/det;
    
    /* dynaK = dynaM*inv(dynaN) */
    dynaK[0][0] = dynaM[0][0]*inv[0][0]+dynaM[1][0]*inv[0][1];
    dynaK[1][0] = dynaM[0][0]*inv[1][0]+dynaM[1][0]*inv[1][1];
    dynaK[0][1] = dynaM[0][1]*inv[0][0]+dynaM[1][1]*inv[0][1];
    dynaK[1][1] = dynaM[0][1]*inv[1][0]+dynaM[1][1]*inv[1][1];
}

void dray_central (float** traj, int length, float* o, float* d, int* n, int* m)
/*< project traced ray onto grid >*/
{
    int it, i, j;
    int tempz, tempy;
    
    for (it=0; it < length; it++) {
	/* NOTE: debug */

	tempz = (int)floorf((traj[it][0]-o[0])/d[0]);
	tempy = (int)floorf((traj[it][1]-o[1])/d[1]);

	for (j=tempy; j <= (tempy+1<n[1]? tempy+1: tempy); j++) {
	    for (i=tempz; i <= (tempz+1<n[0]? tempz+1: tempz); i++) {
		m[j*n[0]+i] = 1;
	    }
	}
/*
	tempz = (int)floorf((traj[it][0]-o[0])/d[0]+0.5);
	tempy = (int)floorf((traj[it][1]-o[1])/d[1]+0.5);
*/	
	m[tempy*n[0]+tempz] = 1;
    }
}

int dray_search (float** traj, int length, float* x)
/*< search for nearest point, return index on central ray >*/
{
    int it, x0;
    double dist, min;

    min = 1.e10;
    x0 = 0;

    /* NOTE: naive search - loop over all possibilities */
    for (it=0; it < length; it++) {
	dist = hypot((double)traj[it][0]-x[0],(double)traj[it][1]-x[1]);

	if (dist < min) {
	    min = dist;
	    x0 = it;
	}
    }

    return x0;
}
