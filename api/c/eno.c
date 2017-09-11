/* 1-D ENO interpolation */
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
#include <float.h>
#include <stdlib.h>

#include "eno.h"
#include "alloc.h"

#include "_defs.h"
/*^*/

#ifndef _sf_eno_h

typedef struct Eno *sf_eno;
/* abstract data type */
/*^*/

typedef enum {FUNC, DER, BOTH} der;
/* flag values */
/*^*/

#endif

struct Eno {
    int order, n;
    float** diff;
};
/* concrete data type */

sf_eno sf_eno_init (int order /* interpolation order */, 
		    int n     /* data size */)
/*< Initialize interpolation object. >*/
{
    sf_eno ent;
    int i;
    
    ent = (sf_eno) sf_alloc(1,sizeof(*ent));
    ent->order = order;
    ent->n = n;
    ent->diff = (float**) sf_alloc(order,sizeof(float*));
    for (i = 0; i < order; i++) {
	ent->diff[i] = sf_floatalloc(n-i);
    }
  
    return ent;
}

void sf_eno_close (sf_eno ent)
/*< Free internal storage >*/
{
    int i;

    for (i = 0; i < ent->order; i++) {
	free(ent->diff[i]);
    }
    free (ent->diff);
    free (ent);
}

void sf_eno_set (sf_eno ent, float* c /* data [n] */)
/*< Set the interpolation table. c can be changed or freed afterwords >*/
{
    int i, j;
    
    for (i=0; i < ent->n; i++) {
	/* copy the initial data */
	ent->diff[0][i] = c[i];
    }
    
    for (j=1; j < ent->order; j++) {
	for (i=0; i < ent->n-j; i++) {
	    /* compute difference tables */
	    ent->diff[j][i] = ent->diff[j-1][i+1] - ent->diff[j-1][i];
	}
    }
}

void sf_eno_set_wstride (sf_eno ent, float* c /* data [n] */, int stride)
/*< Set the interpolation table. c can be changed or freed afterwords >*/
{
    int i, j;

    if (stride < 1) stride = 1;

    for (i=0; i < ent->n; i++) {
	/* copy the initial data */
	ent->diff[0][i] = c[i*stride];
    }
    
    for (j=1; j < ent->order; j++) {
	for (i=0; i < ent->n-j; i++) {
	    /* compute difference tables */
	    ent->diff[j][i] = ent->diff[j-1][i+1] - ent->diff[j-1][i];
	}
    }
}

void sf_eno_apply (sf_eno ent, 
		int i     /* grid location */, 
		float x   /* offset from grid */, 
		float *f  /* output data value */, 
		float *f1 /* output derivative */, 
		der what  /* flag of what to compute */) 
/*< Apply interpolation >*/
{
    int j, k, i1, i2, n;
    float s, s1, y, w, g, g1;
    
    i2 = SF_MAX(0,SF_MIN(i,ent->n-ent->order));
    i1 = SF_MIN(i2,SF_MAX(0,i-ent->order+2));
    
    w = fabsf(ent->diff[ent->order-1][i1]);
    for (j=i1+1; j <=i2; j++) {
	g = fabsf(ent->diff[ent->order-1][j]);
	if (w > g) w = g;
    }
    
    /* loop over starting points */
    for (g = 0., g1 = 0., n = 0, j = i1; j <= i2; j++) {
	if (fabsf(ent->diff[ent->order-1][j]) > w) continue;
	n++;
        
	y = x + i - j;
	
	/* loop to compute the polynomial */
	for (s = 1., s1 = 0., k=0; k < ent->order; k++) {
	    if (what != FUNC) {
		g1 += s1*ent->diff[k][j];
		s1 = (s + s1*(y-k))/(k+1.);
	    }
	    if (what != DER) g += s*ent->diff[k][j];
	    s *= (y-k)/(k+1.);
	}
    }
    
    if (what != DER) *f = g/n;
    if (what != FUNC) *f1 = g1/n;
}

/* 	$Id$	 */

