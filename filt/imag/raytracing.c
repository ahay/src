#include <rsf.h>

#include "raytracing.h"
#include "grid2.h"
#include "grid3.h"
#include "atela.h"

/* concrete data type */
struct RayTrace {
    int dim, nt;
    float dt, z0;
    grid3 grd3;
    grid2 grd2;
};

/*
  Function: raytrace_init
  -----------------------
  Initialize ray tracing object
  dim               - dimensionality (2 or 3)
  nt                - number of ray tracing steps
  dt                - ray tracing step (in time)
  n[dim]            - slowness dimensions
  o[dim], d[dim]    - slowness grid
  slow2[n3*n2*n1]   - slowness squared
  order             - interpolation order
  Note: 
  1. Increasing order increases accuracy but
  decreases efficiency. Recommended values: 3 or 4.
  2. slow2 can be changed or deallocated after
  raytrace_init.
*/
raytrace raytrace_init(int dim, int nt, float dt,
		       int* n, float* o, float* d,
		       float* slow2, int order)
{
    raytrace rt;
    
    rt = (raytrace) sf_alloc (1,sizeof(*rt));
    
    rt->dim = dim;
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

/*
  Function: raytrace_close
  ------------------------
  Free internal storage
*/
void raytrace_close (raytrace rt)
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

/*
  Function: trace_ray
  -------------------
  Trace a ray
  rt               - ray tracing object
  x[dim]           - point location {z,y,x}
  p[dim]           - ray parameter vector
  traj[(nt+1)*dim] - ray trajectory (output)
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
  4. The output code for it = trace_ray(...)
  it=0 - ray traced to the end without leaving the grid
  it>0 - ray exited at the top of the grid
  it<0 - ray exited at the side or bottom of the grid
  The total traveltime along the ray is 
  nt*dt if (it = 0); abs(it)*dt otherwise 
*/
int trace_ray (raytrace rt, float* x, float* p, float** traj)
{
    int it=0;
    
    switch (rt->dim) {
	case 2:
	    it = atela_step (rt->dim, rt->nt, rt->dt, true, x, p, 
			     rt->grd2, 
			     grid2_vgrad, grid2_vel, grid2_term, traj);
	    break;
	case 3:
	    it = atela_step (rt->dim, rt->nt, rt->dt, true, x, p, 
			     rt->grd3, 
			     grid3_vgrad, grid3_vel, grid3_term, traj);
	    break;
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}


