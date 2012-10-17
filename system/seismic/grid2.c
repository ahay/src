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
    float ***vel;
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
    int i, j;
    float v, vgrad[2];
    grid2 grd;

    grd = (grid2) sf_alloc(1,sizeof(*grd));

    grd->n1 = n1; grd->o1 = o1; grd->d1 = d1; 
    grd->n2 = n2; grd->o2 = o2; grd->d2 = d2;

    grd->pnt = sf_eno2_init (order != 0 ? order : 3, n1, n2);
    sf_eno2_set1 (grd->pnt, slow2);
    grd->vel = NULL;
    if (order <= 0) {
        grd->vel = sf_floatalloc3 (3, n1, n2);
        for (j = 0; j < n2; j++) {
            for (i = 0; i < n1; i++) {
                sf_eno2_apply (grd->pnt, i, j, 0., 0., &v, vgrad, BOTH);
                grd->vel[j][i][0] = slow2[j*n1 + i];
                grd->vel[j][i][1] = 0.5*vgrad[0]/d1;
                grd->vel[j][i][2] = 0.5*vgrad[1]/d2;
            }
        }
        sf_eno2_close (grd->pnt);
    }

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

    if (grd->vel) {
        if (i < 0) {
            i = 0; x = 0.0;
        }
        if (i >= (grd->n1 - 1)) {
            i = grd->n1 - 2; x = 1.0;
        }
        if (j < 0) {
            j = 0; y = 0.0;
        }
        if (j >= (grd->n2 - 1)) {
            j = grd->n2 - 2; y = 1.0;
        }
        f = grd->vel[j][i][0]*(1.0 - y)*(1.0 - x) +
            grd->vel[j][i + 1][0]*(1.0 - y)*x +
            grd->vel[j + 1][i][0]*y*(1.0 - x) +
            grd->vel[j + 1][i + 1][0]*y*x;
    } else
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

    if (grd->vel) {
        if (i < 0) {
            i = 0; x = 0.0;
        }
        if (i >= (grd->n1 - 1)) {
            i = grd->n1 - 2; x = 1.0;
        }
        if (j < 0) {
            j = 0; y = 0.0;
        }
        if (j >= (grd->n2 - 1)) {
            j = grd->n2 - 2; y = 1.0;
        }
        grad[0] = grd->vel[j][i][1]*(1.0 - y)*(1.0 - x) +
                  grd->vel[j][i + 1][1]*(1.0 - y)*x +
                  grd->vel[j + 1][i][1]*y*(1.0 - x) +
                  grd->vel[j + 1][i + 1][1]*y*x;
        grad[1] = grd->vel[j][i][2]*(1.0 - y)*(1.0 - x) +
                  grd->vel[j][i + 1][2]*(1.0 - y)*x +
                  grd->vel[j + 1][i][2]*y*(1.0 - x) +
                  grd->vel[j + 1][i + 1][2]*y*x;
    } else {
        sf_eno2_apply(grd->pnt, i, j, x, y, &f, f1, DER);
        grad[0] = 0.5*f1[0]/grd->d1;
        grad[1] = 0.5*f1[1]/grd->d2;
    }
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
    if (grd->vel) {
        free (grd->vel[0][0]); free (grd->vel[0]); free (grd->vel);
    } else
        sf_eno2_close (grd->pnt);
    free (grd);
}

/* 	$Id$	 */
