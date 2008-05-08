/* 1-D Weighted Power3 ENO5 interpolation */
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

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "wpeno.h"
#include "alloc.h"

#include "_defs.h"
/*^*/

#ifndef _sf_wpeno_h

typedef struct Weno *sf_eno;
/* abstract data type */
/*^*/

typedef enum {FUNC, DER, BOTH} der;
/* flag values */
/*^*/

#endif

#define NWORDER 4
#define SICONST 13/14

struct Weno {
    int n;
    float** diff;
};
/* concrete data type */

sf_eno sf_wpeno5_init (int n     /* data size */)
/*< Initialize interpolation object. >*/
{
    sf_eno ent;
    int i;
    
    ent = (sf_eno) sf_alloc(1,sizeof(*ent));
    ent->n = n;
    ent->diff = (float**) sf_alloc(NWORDER,sizeof(float*));
    for (i = 0; i < NWORDER-1; i++) {
	ent->diff[i] = sf_floatalloc(n-i);
    }
    ent->diff[NWORDER-1] = sf_floatalloc(n-NWORDER+1);
  
    return ent;
}

void sf_wpeno5_close (sf_eno ent)
/*< Free internal storage >*/
{
    int i;

    for (i = 0; i < NWORDER; i++) {
	free(ent->diff[i]);
    }
    free (ent->diff);
    free (ent);
}

float pweno3 (float x, float y) 
/*< Powereno3 limiter. >*/
{
    float a, b, mins, pw;

    if (fabsf(y) >= fabsf(x)) {
	a = fabsf(x);
	b = fabsf(y);
	mins = SF_SIG(x);
    } else {
        b = fabsf(x);
	a = fabsf(y);
	mins = SF_SIG(y);
    }
    pw = mins * a * (pow(a,2) + 3.0*pow(b,2)) / pow((a+b),2) );
    return pw;
}

void sf_eno_set (sf_eno ent, float* c /* data [n] */)
/*< Set the interpolation undivided difference table. c can be changed or freed afterwards >*/
{
    int i, j;
    
    for (i=0; i < ent->n; i++) {
	/* copy the initial data */
	ent->diff[0][i] = c[i];
    }
    
    for (j=1; j < NWORDER-1; j++) {
	for (i=0; i < ent->n-j; i++) {
	    /* compute difference tables */
	    ent->diff[j][i] = ent->diff[j-1][i+1] - ent->diff[j-1][i];
	}
    }
    /* Last column is powereno3 limiter of undivided difference */
    ent->diff[NWORDER-1][i] = pweno3(ent->diff[NWORDER-2][i+1],ent->diff[NWORDER-2][i]);
}

void sf_eno_apply (sf_eno ent, 
		int i     /* grid location */, 
		float x   /* offset from grid */, 
		float *f  /* output data value */, 
		float *f1 /* output derivative */, 
		der what  /* flag of what to compute */) 
/*< Apply interpolation >*/
{
    int j, k, i1, i2, n, order;
    float s, s1, y, w, g, g1;
    float a0, a1, a2, eps, si0, si1, si2;
    
    order = NWORDER-1;
    i2 = SF_MAX(0,SF_MIN(i,ent->n-order));
    i1 = SF_MIN(i2,SF_MAX(0,i-order+2));
    
    w = fabsf(ent->diff[NWORDER-1][i1]);
    for (j=i1+1; j <=i2; j++) {
	g = fabsf(ent->diff[NWORDER-1][j]);
	if (w > g) w = g;
    }
    
    /* loop over starting points */
    for (g = 0., g1 = 0., n = 0, j = i1; j <= i2; j++) {
	if (fabsf(ent->diff[NWORDER-1][j]) > w) continue;
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

        /* Smoothness indicators */
	si0 = SICONST;
	si1 = SICONST;
	si2 = SICONST;

        /* Weights */
	a0 = 0.2/(eps+si0);
	a1 = 0.2/(eps+si1);
	a2 = 0.6/(eps+si2);

	/* loop to create convex combination of polynomials */
    }
    
    if (what != DER) *f = g/n;
    if (what != FUNC) *f1 = g1/n;
}


