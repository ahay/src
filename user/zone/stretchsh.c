/* Inverse shifted-linear interpolation */
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

#include "stretchsh.h"

#ifndef _stretchsh_h

typedef struct Map4 *map4;
/* abstract data type */
/*^*/

#endif

struct Map4 {
    int nt, nd, ib, ie;
    float t0,dt, eps;
    int *x; 
    bool *m;
    float **w, *diag, *offd[1];
    sf_bands slv;
};


static const float shifted[] = {0.21};
static const float sh0=0.79;

void sf_shlin_post (int n            /* total trace length */, 
		      int n1           /* start point */, 
		      int n2           /* end point */, 
		      const float* inp /* spline coefficients */, 
		      float* out       /* function values */)
/*< shifted-linear interpolation post-filtering >*/
{
    int i;
    float o2;

    o2 = shifted[0];

    for (i=n1; i < n2; i++) {
	if (n-1==i) {
	    out[i-n1] = sh0*inp[i] + o2*inp[i-1];
	} else {
	    out[i-n1] = sh0*inp[i] + o2*inp[i+1];
	}
    }
}

static void spline2_int (float x, float* w)
/*Linear Interpolation*/
{
    w[0] = 1. - x; 
    w[1] = x; 
}

map4 stretchsh_init (int n1, float o1, float d1 /* regular axis */, 
		    int nd                     /* data samples */, 
		    float eps                  /* regularization */)
/*< initialize >*/
{
    map4 str;
    
    str = (map4) sf_alloc (1, sizeof(*str));

    str->nt = n1; 
    str->t0 = o1; 
    str->dt = d1; 
    str->nd = nd; 
    str->eps = eps;
    
    str->x = sf_intalloc (nd);
    str->m = sf_boolalloc (nd);
    str->w = sf_floatalloc2 (2,nd);
    str->diag = sf_floatalloc (str->nt);
    str->offd[0] = sf_floatalloc (str->nt-1); // Tridiagonal
    str->slv = sf_banded_init (str->nt,1); // Tridiagonal
/*    str->tslv = sf_spline4_init(str->nt);*/

    return str;
}

void stretchsh_define (map4 str, const float* coord /* [nd] */)
/*< set coordinates >*/
{
    int id, ix, i1, i, n1, i2;
    float rx, d, o, *w;
    
    n1 = str->nt;

    d = str->eps*2.; // regularization filter for linear interpolation
    o = -str->eps;

    for (i1 = 0; i1 < n1; i1++) {
	/* regularization */
	str->diag[i1] = d;
	if (i1 < n1-1) str->offd[0][i1] = o;
    }
    
    for (id = 0; id < str->nd; id++) {
	rx = (coord[id] - str->t0)/str->dt - 0.21; // with optimal shift 
	ix = floorf(rx);
	rx -= ix;

	if (ix <= -2 || ix >= n1) {
	    str->m[id] = true; 
	    continue;
	}

	str->x[id] = ix; 
	str->m[id] = false; 
	w = str->w[id];

	spline2_int(rx,w);
	
	i1 = SF_MAX(0,-ix);
	i2 = SF_MIN(2,n1-ix);

	for (i = i1; i < i2; i++) { 
	    str->diag[ix+i] += w[i] * w[i];
	    str->offd[0][ix+i] += w[i] * w[i+1];
	}
    }

    sf_banded_define (str->slv, str->diag, str->offd);
    
    str->ib = -1;
    for (i1 = 0; i1 < n1; i1++) {
	if (str->diag[i1] != d) {
	    str->ib = i1-1; 
	    break;
	}
    }
    
    str->ie = n1+1;
    for (i1 = n1-1; i1 >= 0; i1--) {
	if (str->diag[i1] != d) {
	    str->ie = i1+2;
	    break;
	}
    }
}


void stretchsh_apply (map4 str, 
		     const float* ord /* [nd] */, 
		     float* mod       /* [n1] */)
/*< transform ordinates to model >*/
{
    int id, it, nt,i, i1, i2;
    float *w, *mm;
    
    mm = str->diag;
    nt = str->nt;

    for (it = 0; it < nt; it++) {
	mm[it] = 0.;
    }
    
    for (id = 0; id < str->nd; id++) {
	if (str->m[id]) continue;
	
	it = str->x[id]; 
	w = str->w[id]; 
	
	i1 = SF_MAX(0,-it);
	i2 = SF_MIN(2,nt-it);

	for (i=i1; i < i2; i++) {
	    mm[it+i] += w[i]*ord[id];
	}
    }    

    sf_banded_solve (str->slv, mm);

    for (it = 0; it <= str->ib; it++) {
	mm[it] = 0.;
    }
    
    for (it = str->ie; it < nt; it++) {
	mm[it] = 0.;
    }

    sf_shlin_post(nt,0,nt,mm,mod);

    for (it = 0; it <= str->ib; it++) {
	mod[it] = 0.;
    }
    
    for (it = str->ie; it < nt; it++) {
	mod[it] = 0.;
    }
}

void stretchsh_close (map4 str)
/*< free allocated storage >*/
{
    free (str->x);
    free (str->m);
    free (str->w[0]);
    free (str->w);
    free (str->diag);
    free (str->offd[0]);
    
    sf_banded_close (str->slv);
    free (str);
}


/* 	$Id: stretch4.c 11336 2013-11-13 21:21:13Z sfomel $	 */
