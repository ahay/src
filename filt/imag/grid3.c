#include <math.h>

#include <rsf.h>

#include "grid3.h"
#include "eno3.h"

/* concrete data type */
struct Grid3 {
    eno3 pnt;
    int n1, n2, n3;
    float o1, d1, o2, d2, o3, d3;
};

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
grid3 grid3_init (int n1, float o1, float d1, 
		  int n2, float o2, float d2,
		  int n3, float o3, float d3,
		  float *slow2, int order)
{
    grid3 grd;
    
    grd = (grid3) sf_alloc(1,sizeof(*grd));
    
    grd->n1 = n1; grd->o1 = o1; grd->d1 = d1; 
    grd->n2 = n2; grd->o2 = o2; grd->d2 = d2;
    grd->n3 = n3; grd->o3 = o3; grd->d3 = d3;
    
    grd->pnt = eno3_init (order, n1, n2, n3);
    eno3_set1 (grd->pnt, slow2);
    
    return grd;
}

/*
  Function: grid3_vel
  -------------------
  Extract a value from the grid
  xy[3] - data coordinates
*/
float grid3_vel(void* par, float* xy)
{
    grid3 grd;
    float x, y, z, f, f1[3];
    int i, j, k;
    
    grd = (grid3) par;
    x = (xy[0]-grd->o1)/grd->d1; i = floor(x); x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = floor(y); y -= j;
    z = (xy[2]-grd->o3)/grd->d3; k = floor(z); z -= k;
    
    eno3_apply(grd->pnt, i, j, k, x, y, z, &f, f1, FUNC);
    return f;
}

/*
  Function: grid3_vgrad
  ---------------------
  Extract (1/2 of) gradient values from the grid
  xy[3]   - data coordinates
  grad[3] - gradient (output)
*/
void grid3_vgrad(void* par, float* xy, float* grad)
{
    grid3 grd;
    float x, y, z, f, f1[3];
    int i, j, k;
    
    grd = (grid3) par;
    x = (xy[0]-grd->o1)/grd->d1; i = floor(x); x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = floor(y); y -= j;
    z = (xy[2]-grd->o3)/grd->d3; k = floor(z); z -= k;
    
    eno3_apply(grd->pnt, i, j, k, x, y, z, &f, f1, DER);
    
    grad[0] = 0.5*f1[0]/grd->d1;
    grad[1] = 0.5*f1[1]/grd->d2;
    grad[2] = 0.5*f1[2]/grd->d3;
}

/* 
   Function: grid3_term
   --------------------
   Termination criterion
   returns 0 if xy (data coordinates)
   are inside the grid
*/
int grid3_term (void *par, float* xy)
{
    grid3 grd;
    
    grd = (grid3) par;
    return (xy[0] < grd->o1 || xy[0] > grd->o1 + (grd->n1-1)*grd->d1 || 
	    xy[1] < grd->o2 || xy[1] > grd->o2 + (grd->n2-1)*grd->d2 ||
	    xy[2] < grd->o3 || xy[2] > grd->o3 + (grd->n3-1)*grd->d3);
}

/* 
   Function: grid3_close
   ---------------------
   Free internal storage
*/
void grid3_close(grid3 grd)
{
    eno3_close (grd->pnt);
    free (grd);
}

/* 	$Id: grid3.c,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */

