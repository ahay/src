/* 2-D velocity grid for ray tracing. */
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

#include "grid2.h"

#ifndef _grid2_h

typedef struct Grid2* grid2;
/* abstract data type */
/*^*/

#endif

struct Grid2 {
    sf_eno2 pnt;
    int n1, n2;
    float o1, d1, o2, d2;
};
/* concrete data type */

grid2 grid2_init (int n1, float o1, float d1 /* first axis */, 
		  int n2, float o2, float d2 /* second axis */,
		  float *slow2               /* data values [n1*n2] */, 
		  int order                  /* interpolation order */)
/*< Initialize grid object >*/
{
    grid2 grd;
    
    grd = (grid2) sf_alloc(1,sizeof(*grd));

    grd->n1 = n1; grd->o1 = o1; grd->d1 = d1; 
    grd->n2 = n2; grd->o2 = o2; grd->d2 = d2;
    
    grd->pnt = sf_eno2_init (order, n1, n2);
    sf_eno2_set1 (grd->pnt, slow2);

    return grd;
}

float grid2_vel(void* par /* grid */, 
		float* xy /* coordinate [2] */)
/*<  Extract a value from the grid >*/
{
    grid2 grd;
    float x, y, f, f1[2];
    int i, j;
    
    grd = (grid2) par;
    x = (xy[0]-grd->o1)/grd->d1; i = x; x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = y; y -= j;
    
    sf_eno2_apply(grd->pnt, i, j, x, y, &f, f1, FUNC);
    return f;
}

void grid2_vgrad(void* par   /* grid */, 
		 float* xy   /* coordinate [2] */, 
		 float* grad /* output gradient [2] */)
/*< Extract (1/2 of) gradient values from the grid >*/
{
    grid2 grd;
    float x, y, f, f1[2];
    int i, j;
    
    grd = (grid2) par;
    x = (xy[0]-grd->o1)/grd->d1; i = floor(x); x -= i;
    y = (xy[1]-grd->o2)/grd->d2; j = floor(y); y -= j;
    
    sf_eno2_apply(grd->pnt, i, j, x, y, &f, f1, DER);
    
    grad[0] = 0.5*f1[0]/grd->d1;
    grad[1] = 0.5*f1[1]/grd->d2;
}

int grid2_term (void* par /* grid */, 
		float* xy /* location [2] */)
/*< Termination criterion. returns 0 if xy (data coordinates)
  are inside the grid >*/
{
    grid2 grd;
    
    grd = (grid2) par;
    return (xy[0] < grd->o1 || xy[0] > grd->o1 + (grd->n1-1)*grd->d1 || 
	    xy[1] < grd->o2 || xy[1] > grd->o2 + (grd->n2-1)*grd->d2);
}

void grid2_close(grid2 grd)
/*< Free internal storage >*/
{
    sf_eno2_close (grd->pnt);
    free (grd);
}

/* 	$Id: grid2.c 7107 2011-04-10 02:04:14Z ivlad $	 */
