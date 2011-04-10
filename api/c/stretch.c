/* Inverse linear interpolation */
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

#include "stretch.h"

#include "_bool.h"
/*^*/

#include "tridiagonal.h"
#include "alloc.h"

#ifndef _sf_stretch_h

typedef struct sf_Map *sf_map;
/* abstract data type */
/*^*/

#endif


struct sf_Map {
    int nt, nd, ib, ie;
    float t0,dt, eps;
    int *x; 
    bool *m, narrow;
    float *w, *diag, *offd;
    sf_tris slv;
};

sf_map sf_stretch_init (int n1, float o1, float d1 /* regular axis */, 
			int nd                     /* data length */, 
			float eps                  /* regularization */, 
			bool narrow                /* if zero boundary */)
/*< initialize >*/
{
    sf_map str;
    
    str = (sf_map) sf_alloc (1, sizeof(*str));

    str->nt = n1; 
    str->t0 = o1; 
    str->dt = d1; 
    str->nd = nd; 
    str->eps = eps;
    str->narrow = narrow;
    
    str->x = sf_intalloc (nd);
    str->m = sf_boolalloc (nd);
    str->w = sf_floatalloc (nd);
    str->diag = sf_floatalloc (n1);
    str->offd = sf_floatalloc (n1-1);
  
    str->slv = sf_tridiagonal_init (n1);

    return str;
}

void sf_stretch_define (sf_map str, const float* coord)
/*< define coordinates >*/
{
    int id, ix, i1;
    float rx, r1;
    
    for (i1 = 0; i1 < str->nt-1; i1++) {
	/* regularization */
	str->diag[i1] = str->eps;
	str->offd[i1] = -0.5*str->eps;
    }
    if (str->narrow) {
	str->diag[str->nt-1] = str->eps;
    } else {
	str->diag[0] = 0.5*str->eps;
	str->diag[str->nt-1] = 0.5*str->eps;
    }
    
    for (id = 0; id < str->nd; id++) {
	rx = (coord[id] - str->t0)/str->dt; 
	ix = floorf(rx); 
	rx = rx - ix;
	if (ix < 0 || ix > str->nt - 2) {
	    str->m[id] = true; 
	    continue;
	}
	str->x[id] = ix; 
	str->m[id] = false; 
	str->w[id] = rx;
	
	r1 = 1. - rx;
	
	str->diag[ix]   += r1 * r1;
	str->diag[ix+1] += rx * rx;
	str->offd[ix]   += r1 * rx;
    }
    
    sf_tridiagonal_define (str->slv, str->diag, str->offd);
    
    if (str->narrow) {
	str->ib = -1;
	for (i1 = 0; i1 < str->nt; i1++) {
	    if (str->diag[i1] != str->eps) {
		str->ib = i1-1; 
		break;
	    }
	}

	str->ie = str->nt;
	for (i1 = str->nt-1; i1 >= 0; i1--) {
	    if (str->diag[i1] != str->eps) {
		str->ie = i1+1;
		break;
	    }
	}
    }
}

void sf_stretch_apply (sf_map str, const float* ord, float* mod)
/*< convert ordinates to model >*/
{
    int id, i1, i2;
    float w1, w2;
    
    for (i1 = 0; i1 < str->nt; i1++) {
	mod[i1] = 0.;
    }
    
    for (id = 0; id < str->nd; id++) {
	if (str->m[id]) continue;
	
	i1 = str->x[id]; i2 = i1+1;
	w2 = str->w[id]; w1 = 1.-w2;
	
	mod[i1] += w1 * ord[id];
	mod[i2] += w2 * ord[id];
    }
    
    sf_tridiagonal_solve (str->slv, mod);

    if (str->narrow) {
	for (i1 = 0; i1 <= str->ib; i1++) {
	    mod[i1] = 0.;
	}
      
	for (i1 = str->ie; i1 < str->nt; i1++) {
	    mod[i1] = 0.;
	}
    }
}

void sf_stretch_invert (sf_map str, float* ord, const float* mod)
/*< convert model to ordinates by linear interpolation >*/
{
    int id, i1, i2;
    float w1, w2;
    
    for (id = 0; id < str->nd; id++) {
	if (str->m[id]) continue;
	
	i1 = str->x[id]; i2 = i1+1;
	w2 = str->w[id]; w1 = 1.-w2;
	
	ord[id] = w1*mod[i1] + w2*mod[i2];
    }
}

void sf_stretch_close (sf_map str)
/*< free allocated storage >*/
{
    free (str->x);
    free (str->m);
    free (str->w);
    free (str->diag);
    free (str->offd);
    
    sf_tridiagonal_close (str->slv);
    free (str);
}

/* 	$Id$	 */
