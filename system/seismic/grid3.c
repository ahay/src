/* 3-D velocity grid for ray tracing. */
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

#include <math.h>

#include <rsf.h>

#include "grid3.h"
#include "eno3.h"

#ifndef _grid3_h

typedef struct Grid3* grid3;
/* abstract data type */
/*^*/

#endif

struct Grid3 {
    eno3 pnt;
    int n1, n2, n3;
    float o1, d1, o2, d2, o3, d3;
};
/* concrete data type */

grid3 grid3_init (int n1, float o1, float d1 /* first axis */, 
		  int n2, float o2, float d2 /* second axis */,
		  int n3, float o3, float d3 /* third axis */,
		  float *slow2               /* data [n1*n2*n3] */, 
		  int order                  /* interpolation order */)
/*< Initialize 3-D grid. >*/
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

float grid3_vel(void* par /* grid */, 
		float* xy /* location [3] */)
/*< Extract a value from the grid. >*/
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

void grid3_vgrad(void* par /* grid */, 
		 float* xy /* location [3] */, 
		 float* grad /* output gradient [3] */)
/*< Extract (1/2 of) gradient values from the grid >*/
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

int grid3_term (void *par /* grid */, 
		float* xy /* location [3] */)
/*< Termination criterion. Returns 0 if xy (data coordinates)
  are inside the grid >*/
{
    grid3 grd;
    
    grd = (grid3) par;
    return (xy[0] < grd->o1 || xy[0] > grd->o1 + (grd->n1-1)*grd->d1 || 
	    xy[1] < grd->o2 || xy[1] > grd->o2 + (grd->n2-1)*grd->d2 ||
	    xy[2] < grd->o3 || xy[2] > grd->o3 + (grd->n3-1)*grd->d3);
}

void grid3_close(grid3 grd)
/*< Free internal storage >*/
{
    eno3_close (grd->pnt);
    free (grd);
}

/* 	$Id$	 */

