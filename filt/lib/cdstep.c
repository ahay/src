/* Conjugate-direction iteration */
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

#include <float.h>

#include "cdstep.h"
#include "llist.h"
#include "alloc.h"

#include "_bool.h"
/*^*/

static sf_list steps;

static double dotprod (int n, const float* x, const float* y); 
/* double-precision dot product */

void sf_cdstep_init(void) 
/*< initialize internal storage >*/
{
    steps = sf_llist_init();
}

void sf_cdstep_close(void) 
/*< free internal storage >*/
{
    sf_llist_close(steps);
}

void sf_cdstep(bool forget     /* restart flag */, 
	       int nx          /* model size */, 
	       int ny          /* data size */, 
	       float* x        /* current model [nx] */, 
	       const float* g  /* gradient [nx] */, 
	       float* rr       /* data residual [ny] */, 
	       const float* gg /* conjugate gradient [ny] */) 
/*< Step of conjugate-direction iteration. 
  The data residual is rr = A x - dat
>*/
{
    float *s, *si;
    double alpha, beta;
    int i, n, ix, iy;

    s = sf_floatalloc(nx+ny);

    for (ix=0; ix < nx; ix++) {
	s[ix] = g[ix];
    }
    for (iy=0; iy < ny; iy++) {
	s[nx+iy] = gg[iy];
    }

    sf_llist_rewind (steps);
    n = sf_llist_depth (steps);

    for (i=0; i < n; i++) {
	sf_llist_down (steps, &si, &beta);
	alpha = dotprod (ny, gg, si+nx) / beta;
	for (ix=0; ix < nx+ny; ix++) {
	    s[ix] -= alpha * si[ix];
	}
    }

    beta = dotprod (ny, s+nx, s+nx);
    if (beta < DBL_EPSILON) return;

    sf_llist_add (steps, s, beta);
    if (forget) {
	sf_llist_chop (steps);
    }
    alpha = - dotprod (ny, rr, s+nx) / beta;
    
    for (ix=0; ix < nx; ix++) {
	x[ix] += alpha * s[ix];
    }
    for (iy=0; iy < ny; iy++) {
	rr[iy] += alpha * s[nx+iy];
    }
}

void sf_cdstep_diag(int nx, float *res /* [nx] */)
/*< compute diagonal of the model resolution matrix >*/
{
    int i, n, ix;
    float *s;
    double sn;

    sf_llist_rewind (steps);
    n = sf_llist_depth (steps);

    for (ix=0; ix < nx; ix++) {
	res[ix] = 0.;
    }

    for (i=0; i < n; i++) {
	sf_llist_down (steps, &s, &sn);
	for (ix=0; ix < nx; ix++) {
	    res[ix] += dotprod(nx,s,s)/sn;
	}
    }
}

void sf_cdstep_mat (int nx, float **res /* [nx][nx] */)
/*< compute complete model resolution matrix >*/
{
    int i, n, i1, i2;
    float *s;
    double sn;

    sf_llist_rewind (steps);
    n = sf_llist_depth (steps);

    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nx; i1++) {
	    res[i2][i1] = 0.;
	}
    }

    for (i=0; i < n; i++) {
	sf_llist_down (steps, &s, &sn);
	for (i2=0; i2 < nx; i2++) {
	    for (i1=0; i1 < nx; i1++) {
		res[i2][i1] += s[i2]*s[i1]/sn;
	    }
	}
    }
}

static double dotprod (int n, const float* x, const float* y) 
/* double-precision dot product */
{
    double prod;
    int i;
    prod = 0.;
    for (i = 0; i < n; i++) {
	prod += ((double) x[i])*y[i];
    }
    return prod;
}
