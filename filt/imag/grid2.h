#ifndef _grid2_h
#define _grid2_h

/*
  File: grid2.h
  -------------
  Operations on 2-D regular grids
*/

/* abstract data type */
typedef struct Grid2* grid2;

/*
  Function: grid2_init
  --------------------
  Initialize 2-D grid
  n1, n2        - grid dimensions
  o1, o2        - grid coordinates
  d1, d2        - grid spacing
  slow2[n2][n1] - data values
  order         - interpolation order
*/
grid2 grid2_init (int nx, float x0, float dx, 
		  int ny, float y0, float dy,
		  float *slow2, int order);

/*
  Function: grid2_vel
  -------------------
  Extract a value from the grid
  xy[2] - data coordinates
*/
float grid2_vel(void* par, float* xy);

/*
  Function: grid2_vgrad
  ---------------------
  Extract (1/2 of) gradient values from the grid
  xy[2]   - data coordinates
  grad[2] - gradient (output)
*/
void grid2_vgrad(void* par, float* xy, float* grad);

/* 
   Function: grid2_term
   --------------------
   Termination criterion
   returns 0 if xy (data coordinates)
   are inside the grid
*/
int grid2_term(void* par, float* xy);

/* 
   Function: grid2_close
   ---------------------
   Free internal storage
*/
void grid2_close(grid2 grd);

#endif

/* 	$Id: grid2.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
