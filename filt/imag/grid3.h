#ifndef _grid3_h
#define _grid3_h

/*
  File: grid3.h
  -------------
  Operations on 3-D regular grids
*/

/* abstract data type */
typedef struct Grid3* grid3;

/*
  Function: grid3_init
  --------------------
  Initialize 3-D grid
  n1, n2, n3        - grid dimensions
  o1, o2, o3        - grid coordinates
  d1, d2, d3        - grid spacing
  slow2[n3][n2][n1] - data values
  order             - interpolation order
*/
grid3 grid3_init (int nz, float z0, float dz, 
		  int ny, float y0, float dy,
		  int nx, float x0, float dx,
		  float *slow2, int order);

/*
  Function: grid3_vel
  -------------------
  Extract a value from the grid
  xy[3] - data coordinates
*/
float grid3_vel(void* par, float* x);

/*
  Function: grid3_vgrad
  ---------------------
  Extract (1/2 of) gradient values from the grid
  xy[3]   - data coordinates
  grad[3] - gradient (output)
*/
void grid3_vgrad(void* par, float* x, float* grad);

/* 
   Function: grid3_term
   --------------------
   Termination criterion
   returns 0 if xy (data coordinates)
   are inside the grid
*/
int grid3_term(void* par, float* x);

/* 
   Function: grid3_close
   ---------------------
   Free internal storage
*/
void grid3_close(grid3 grd);

#endif

/* 	$Id: grid3.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
