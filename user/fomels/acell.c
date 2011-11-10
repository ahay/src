/* Adaptive gridding. */
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

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <rsf.h>

#include "acell.h"

#ifndef _acell_h

typedef struct ACell *acell;
/* abstract data type */
/*^*/

#endif

struct ACell {
    int parent, level;
    float offset;
    acell child[2];
    float* node[2];
};

static int nd, maxsplit, itmp;
static float **tmp;

static void interp_lin (acell cell, float x, float* f);
/* linear interpolation */

static void free_children (acell cell);

static void untree (acell cell, acell *flat);
/* flatten the grid */

acell acell_init (void) 
/*< initialize an adaptive cell >*/
{
    acell cell;
    
    cell = (acell) sf_alloc (1,sizeof(*cell));

    cell->parent = 0;
    cell->offset = 0.;
    cell->level = 1;

    return cell;
}

void acell_set (acell cell, float* dat1 /* left */, float* dat2 /* right */) 
/*< set the left and right points in a cell >*/
{
    cell->node[0] = dat1;
    cell->node[1] = dat2;
}

void free_acell (acell cell) 
/*< free allocated storage >*/
{
    free_children (cell);
    free (cell);
}

void acells_init (int n   /* number of points */, 
		  int max /* maximum splitting */) 
/*< initialize a grid of adaptive cells >*/
{
    nd=n;
    maxsplit=max;
    tmp = sf_floatalloc2(nd,maxsplit);
}

void free_acells (void)
/*< free adaptive cells >*/
{
  free (tmp[0]);
  free (tmp);
}

static void interp_lin (acell cell, float x, float* f) 
/* linear interpolation */
{
    float f1, f2;
    int i;

    assert (x <= 1. && x >= 0.);

    for (i=0; i < nd; i++) {
	f1 = cell->node[0][i];
	f2 = cell->node[1][i] - f1;
	f[i] = f1 + x*f2;
    }
}

void interp_acell (acell cell, float x, float* f) 
/*< interpolate from adaptive cell at position x >*/
{
    int p, q, k, i;
    acell c;

    interp_lin (cell, x, tmp[0]);

    if (x > 0. && x < 1.) {
	for (p=0, c = cell; c->parent > 0; p++) {
	    if (x < 0.5) {
		x = 2*x;
		c = c->child[0];
	    } else {
		x = 2*x-1;
		c = c->child[1];
	    }
	    interp_lin (c, x, tmp[p+1]);
/*	    fprintf(stderr,"%d? %f %f %f\n",p,tmp[p+1][0],tmp[p+1][1]); */
	    if (x > 0. && x < 1.) {
		for (q=p, k=4; q >= 0; q--, k *= 2) {
		    for (i=0; i < nd; i++) {
			tmp[q][i] = tmp[q+1][i]+(tmp[q+1][i]-tmp[q][i])/(k-1.);
		    }
		}
/*		fprintf(stderr,"%d: %f %f\n",p,tmp[0][0],tmp[0][1]); */
	    } else {
		memcpy (f,tmp[p+1],nd*sizeof(float));
		return;
	    }
	}
    }

    memcpy (f,tmp[0],nd*sizeof(float));
}

void fill_cell (acell cell, int i, void* dat, 
		void (*fill)(float x, void* dat, float* f)) 
/*< fill cell and children using fill function >*/
{
    if (cell->parent == 0) {
	fill(i+cell->offset,dat,cell->node[0]);
    } else {
	fill_cell (cell->child[0], i, dat, fill);
	fill_cell (cell->child[1], i, dat, fill);
    }
}

static void free_children (acell cell) 
{
    if (cell->parent > 0) {
	free (cell->child[0]->node[1]);
	free_children(cell->child[0]);
	free_children(cell->child[1]);
	free (cell->child[0]);
	free (cell->child[1]);
	cell->parent = 0;
    }
}

void split_cell (acell cell, 
		 float min1, float max1, 
		 float min2, float max2, int i, void* dat,
		 void (*fill)(float x, void* dat, float* f)) 
/*< split adaptive cell >*/
{
    float dx, dy, dp, dq;

    if (cell->parent>0) {
	dx = fabs(cell->node[0][0]   -cell->child[0]->node[1][0]);
	dy = fabs(cell->node[1][0]   -cell->child[0]->node[1][0]);
	dp = fabs(cell->node[0][nd-1]-cell->child[0]->node[1][nd-1]);
	dq = fabs(cell->node[1][nd-1]-cell->child[0]->node[1][nd-1]);
	if (dx < min1 && dy < min1 && dp < min2 && dq < min2) {
	    free_children (cell);
	} else {
	    split_cell (cell->child[0], min1, max1, min2, max2, i, dat, fill);
	    split_cell (cell->child[1], min1, max1, min2, max2, i, dat, fill);
	}
    } else if (cell->level < maxsplit) {
	dx = fabs(cell->node[0][0]   -cell->node[1][0]);
	dp = fabs(cell->node[0][nd-1]-cell->node[1][nd-1]);
	if (dx > max1 || dp > max2) { /* split */
	    cell->parent = 1;
	    cell->child[0] = (acell) sf_alloc (1,sizeof(struct ACell));
	    cell->child[1] = (acell) sf_alloc (1,sizeof(struct ACell));
	    cell->child[0]->parent = 0;
	    cell->child[1]->parent = 0;
	    cell->child[0]->level = cell->level+1;
	    cell->child[1]->level = cell->level+1;
	    cell->child[0]->offset = cell->offset;
	    cell->child[1]->offset = cell->offset+powf(0.5,cell->level);
	    cell->child[0]->node[0] = cell->node[0];
	    cell->child[0]->node[1] = sf_floatalloc (nd);
	    cell->child[1]->node[0] = cell->child[0]->node[1];
	    cell->child[1]->node[1] = cell->node[1];
	    fill (i+cell->child[1]->offset,dat,cell->child[1]->node[0]);
	    split_cell (cell->child[0], min1, max1, min2, max2, i, dat, fill);
	    split_cell (cell->child[1], min1, max1, min2, max2, i, dat, fill);
	}
    } 
}

void fill_node (acell cell, int i, float x, void* dat,
		void (*fill)(float x, void* dat, float* f)) 
/*< fill node >*/
{
    fill(x,dat,cell->node[i]);
}

float* get_node (acell cell, int i) 
/*< extract node value >*/
{ 
    return cell->node[i];
}

int cell_size (acell cell) 
/*< return cell size >*/
{
    int size;
    if (cell->parent == 0) {
	size = 1;
    } else {
	size = 
	    cell_size (cell->child[0]) +
	    cell_size (cell->child[1]);
    }
    return size;
}


void flat_cell (acell cell, acell *flat, float** dat) 
/*< flatten the grid >*/
{
    int i;

    itmp=0;
    untree (cell,flat);
    for (i=0; i < itmp; i++) {
	dat[i] = flat[i]->node[0];
    }
}

static void untree (acell cell, acell *flat)
/* flatten the cell */
{
    if (cell->parent == 0) {
	flat[itmp] = cell;
	itmp++;
    } else {
	untree (cell->child[0],flat);
	untree (cell->child[1],flat);
    }
}

/* 	$Id: acell.c 854 2004-11-05 16:37:54Z fomels $	 */
