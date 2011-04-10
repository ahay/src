/* 1-D ENO power-p interpolation */
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

#include "pweno.h"
#include "alloc.h"

#include "_defs.h"
/*^*/

#ifndef _sf_pweno_h

typedef struct Pweno *sf_pweno;
/* abstract data type */
/*^*/

typedef enum {FUNC1, DER1, BOTH1} derr;
/* flag values */
/*^*/

#endif

struct Pweno {
    int order, n;
    float** diff;
};
/* concrete data type */

sf_pweno sf_pweno_init (int order /* interpolation order */,
              int n     /* data size */)
/*< Initialize interpolation object >*/
{
    sf_pweno ent;
    int i;
    
    ent = (sf_pweno) sf_alloc(1,sizeof(*ent));
    ent->order = order;
    ent->n = n;
    ent->diff = (float**) sf_alloc(order+1,sizeof(float*));
    for (i = 0; i < order; i++) {
	ent->diff[i] = sf_floatalloc(n-i);
    }
    ent->diff[order] = sf_floatalloc(n-order+1);
  
    return ent;
}

void sf_pweno_close (sf_pweno ent)
/*< Free internal storage >*/
{
    int i;

    for (i = 0; i < ent->order + 1; i++) {
	free(ent->diff[i]);
    }
    free (ent->diff);
    free (ent);
}

float powerpeno (float x, float y, int p /* power order */)
/*< Limiter power-p eno >*/
{
    float s, d, mins, power;

    s = x+y;
    d = x-y;
    mins = ((fabs(y) >= fabs(x)) ? SF_SIG(x) : SF_SIG(y));
    power = 0.0;
    if (fabs(s) > 0.0) power = 0.5*s*(1.0-fabsf(pow(d/s,p)));
    return (mins * power);
}

void sf_pweno_set (sf_pweno ent, float* c /* data [n] */, int p /* power order */)
/*< Set the interpolation undivided difference table. c can be changed or freed afterwards >*/
{
    int i, j;
    
    for (i=0; i < ent->n; i++) {
	/* copy the initial data */
	ent->diff[0][i] = c[i];
    }
    
    for (j=1; j < ent->order; j++) {
	for (i=0; i < ent->n-j; i++) {
	    /* compute undivided difference tables */
	    ent->diff[j][i] = ent->diff[j-1][i+1] - ent->diff[j-1][i];
	}
    }

    for (i=0; i < ent->n-ent->order+1; i++) {
        ent->diff[ent->order][i] = ent->diff[ent->order-1][i];}
    /* Power-p eno limiter of highest order undivided difference */
    for (i=0; i < ent->n-ent->order+1; i++) {
        ent->diff[ent->order-1][i] = powerpeno(ent->diff[ent->order][i+1],ent->diff[ent->order][i],p);}
}

void sf_pweno_apply (sf_pweno ent, 
		int i     /* grid location */, 
		float x   /* offset from grid */, 
		float *f  /* output data value */, 
		float *f1 /* output derivative */, 
		derr what /* flag of what to compute */) 
/*< Apply interpolation >*/
{
    int j, k, i1, i2, n;
    float s, s1, y, w, g, g1, tw;

    /* Stencil */
    i2 = SF_MAX(0,SF_MIN(i,ent->n-ent->order));
    i1 = SF_MIN(i2,SF_MAX(0,i-ent->order+2));
    
    /* ENO selection */
    w = fabsf(ent->diff[ent->order-1][i1]);
    for (j=i1+1; j <=i2; j++) {
	g = ((j != i) ? fabsf(ent->diff[ent->order-1][j]) : fabsf(ent->diff[ent->order][i]) );
	if (w > g) w = g;
    }
    
    /* Loop over starting points */
    for (g = 0., g1 = 0., n = 0, j = i1; j <= i2; j++) {
        tw = ((j != i) ? ent->diff[ent->order-1][j] : ent->diff[ent->order][i] );
	if (fabsf(tw) > w) continue;
	n++;
        
	y = x + i - j;
	
	/* Loop to compute the polynomial */
	for (s = 1., s1 = 0., k=0; k < ent->order; k++) {
	    if (what != FUNC1) {
		g1 += s1*ent->diff[k][j];
		s1 = (s + s1*(y-k))/(k+1.);
	    }
	    if (what != DER1) g += s*((k != (ent->order-1) || j != i) ?  ent->diff[k][j] : ent->diff[ent->order][i] );
	    s *= (y-k)/(k+1.);
	}
    }
    
    if (what != DER1) *f = g/n;
    if (what != FUNC1) *f1 = g1/n;
}


