#include <math.h>

#include <rsf.h>

#include "grid2.h"
#include "eno2.h"

/* concrete data type */
struct Grid2 {
    eno2 pnt;
    int n1, n2;
    float o1, d1, o2, d2;
};

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
grid2 grid2_init (int n1, float o1, float d1, 
		  int n2, float o2, float d2,
		  float *slow2, int order)
{
    grid2 grd;
    
    grd = (grid2) sf_alloc(1,sizeof(*grd));

    grd->n1 = n1; grd->o1 = o1; grd->d1 = d1; 
    grd->n2 = n2; grd->o2 = o2; grd->d2 = d2;
    
    grd->pnt = eno2_init (order, n1, n2);
    eno2_set1 (grd->pnt, slow2);

    return grd;
}

/*
  Function: grid2_vel
  -------------------
  Extract a value from the grid
  xy[2] - data coordinates
*/
float grid2_vel(void* par, float* xy)
{
    grid2 grd;
    float x, y, f, f1[2];
    int i, j;
    
    grd = (grid2) par;
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    
    eno2_apply(grd->pnt, i, j, x, y, &f, f1, FUNC);
    return f;
}

/*
  Function: grid2_vgrad
  ---------------------
  Extract (1/2 of) gradient values from the grid
  xy[2]   - data coordinates
  grad[2] - gradient (output)
*/
void grid2_vgrad(void* par, float* xy, float* grad)
{
    grid2 grd;
    float x, y, f, f1[2];
    int i, j;
    
    grd = (grid2) par;
    x = (xy[0]-grd->o1)/grd->d1; i = floor(x); x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = floor(y); y -= j;
    
    eno2_apply(grd->pnt, i, j, x, y, &f, f1, DER);
    
    grad[0] = 0.5*f1[0]/grd->d1;
    grad[1] = 0.5*f1[1]/grd->d2;
}

/* 
   Function: grid2_term
   --------------------
   Termination criterion
   returns 0 if xy (data coordinates)
   are inside the grid
*/
int grid2_term (void* par, float* xy)
{
    grid2 grd;
    
    grd = (grid2) par;
    return (xy[0] < grd->o1 || xy[0] > grd->o1 + (grd->n1-1)*grd->d1 || 
	    xy[1] < grd->o2 || xy[1] > grd->o2 + (grd->n2-1)*grd->d2);
}

/* 
   Function: grid2_close
   ---------------------
   Free internal storage
*/
void grid2_close(grid2 grd)
{
    eno2_close (grd->pnt);
    free (grd);
}

/* 	$Id: grid2.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
