/* 1-D adaptive grid */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "grid1.h"

#ifndef _grid1_h

typedef struct Grid1 *grid1;
/* abstract data type */
/*^*/

#endif

struct Grid1 {
    float x[2];
    float *val[2];
    grid1 child[2];
};
/* concrete data type */

grid1 grid1_init (void)
/*< constructor >*/
{
    grid1 grid;

    grid = (grid1) sf_alloc (1,sizeof(struct Grid1));
    grid->x[0] = -SF_HUGE;
    grid->x[1] = +SF_HUGE;
    grid->val[0] = (float*) NULL;
    grid->val[1] = (float*) NULL;
    grid->child[0] = (grid1) NULL;
    grid->child[1] = (grid1) NULL;

    return grid;
}

static grid1 grid1_locate(grid1 grid, float coord)
/* locate an interval by binary search */
{
    grid1 left;

    left =  grid->child[0];
    if (NULL == left) return grid;
    if (coord < left->x[1]) return grid1_locate(left,coord);
    return grid1_locate(grid->child[1],coord);
}

void grid1_insert (grid1 grid, float coord, int nv, const float *value)
/*< insert a value in a grid >*/
{
    float *val;
    grid1 node, left, right;

    node = grid1_locate(grid,coord);
    val = sf_floatalloc(nv);
    memcpy(val,value,nv*sizeof(float));

    left = node->child[0] = grid1_init();
    left->x[0] = node->x[0];
    left->x[1] = coord;
    left->val[0] = node->val[0];
    left->val[1] = val;

    right = node->child[1] = grid1_init();
    right->x[0] = coord;
    right->x[1] = node->x[1];
    right->val[0] = val;
    right->val[1] = node->val[1];
}

void grid1_close (grid1 grid)
/*< free allocated storage including values >*/
{
    if (NULL != grid->child[0]) {
	grid1_close(grid->child[0]);
	grid1_close(grid->child[1]);
    } else if (NULL != grid->val[0]) {
	free (grid->val[0]);
	grid->val[0] = NULL;
    }
    free(grid);
}

void grid1_write (grid1 grid, int nv, sf_file out)
/*< write it out >*/
{
    if (NULL != grid->child[0]) {
	grid1_write(grid->child[0],nv,out);
	grid1_write(grid->child[1],nv,out);
    } else if (NULL != grid->val[0]) {
	sf_floatwrite(grid->x,1,out);
	sf_floatwrite(grid->val[0],nv,out);
    }
}

void grid1_interpolate (grid1 grid, float coord, int nv, float *value)
/*< interpolate a value from a grid >*/
{
    int iv;
    float fl, fr, *vl, *vr;
    grid1 node;

    node = grid1_locate(grid,coord);

    vl = node->val[0];
    vr = node->val[1];

    /* linear interpolation -> change later to recursive */
    if (NULL == vl) {
	if (NULL == vr) sf_error ("%s: No values in the grid",__FILE__);

	for (iv=0; iv < nv; iv++) {
	    value[iv] = vr[iv];
	}
    } else if (NULL == vr) {
	if (NULL == vl) sf_error ("%s: No values in the grid",__FILE__);

	for (iv=0; iv < nv; iv++) {
	    value[iv] = vl[iv];
	}
    } else {
	/* linear interpolation */
	fr = (coord - node->x[0])/(node->x[1] - node->x[0]);
	fl = (node->x[1] - coord)/(node->x[1] - node->x[0]);
	/* danger - division by zero */

	for (iv=0; iv < nv; iv++) {
	    value[iv] = fl*vl[iv] + fr*vr[iv];
	}
    }
}
