#ifndef _atela_h
#define _atela_h

#include <rsf.h>

/* 
   File: atela.h
   -------------
   Symplectic Runge-Kutta ray tracing.
*/

/*
  Function: atela_step
  --------------------
  Ray tracing
  dim    - dimensionality
  nt     - number of ray tracing steps
  dt     - ray tracing step 
  intime - if step in time (not sigma)
  x[dim] - point location {z,y,x}
  p[dim] - ray parameter vector
  par    - parameters passed to slowness access functions
  vgrad  - function returning 1/2*(gradient of slowness squared)
  slow2  - function returning slowness squared
  term   - function returning 1 (non-zero) if the ray needs to terminate
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
  it=0 - ray traced to the end without termination
  it>0 - ray terminated
  The total traveltime along the ray is 
  nt*dt if (it = 0); abs(it)*dt otherwise 
*/
int atela_step (int dim, int nt, float dt, bool intime, float* x, float* p, 
		void* par,
		void (*vgrad)(void*,float*,float*), 
		float (*slow2)(void*,float*), 
		int (*term)(void*,float*), float** traj);

#endif
