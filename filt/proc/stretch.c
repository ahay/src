#include <math.h>

#include <rsf.h>

#include "stretch.h"
#include "tridiagonal.h"

struct Map {
    int nt, nd, ib, ie;
    float t0,dt, eps;
    int *x; 
    bool *m, narrow;
    float *w, *diag, *offd;
    tris slv;
};

map stretch_init (int n1, float o1, float d1, int nd, float eps, bool narrow)
{
    map str;
    
    str = (map) sf_alloc (1, sizeof(*str));

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
  
    str->slv = tridiagonal_init (n1);

    return str;
}

void stretch_define (map str, float* coord)
{
    int id, ix, i1;
    float rx, r1;
    
    for (i1 = 0; i1 < str->nt-1; i1++) {
	/* regularization */
	str->diag[i1] = str->eps;
	str->offd[i1] = -0.5*str->eps;
    }
    if (! str->narrow) {
	str->diag[0] = 0.5*str->eps;
	str->diag[str->nt-1] = 0.5*str->eps;
    }
    
    for (id = 0; id < str->nd; id++) {
	rx = (coord[id] - str->t0)/str->dt; 
	ix = floor(rx); 
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
    
    tridiagonal_define (str->slv, str->diag, str->offd);
    
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

void stretch_apply (map str, float* ord, float* mod)
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
    
    tridiagonal_solve (str->slv, mod);

    if (str->narrow) {
	for (i1 = 0; i1 <= str->ib; i1++) {
	    mod[i1] = 0.;
	}
      
	for (i1 = str->ie; i1 < str->nt; i1++) {
	    mod[i1] = 0.;
	}
    }
}

void stretch_invert (map str, float* ord, float* mod)
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

void stretch_close (map str)
{
    free (str->x);
    free (str->m);
    free (str->w);
    free (str->diag);
    free (str->offd);
    
    tridiagonal_close (str->slv);
    free (str);
}

/* 	$Id: stretch.c,v 1.4 2004/04/03 02:41:17 fomels Exp $	 */

