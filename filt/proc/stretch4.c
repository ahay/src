#include <math.h>

#include <rsf.h>

#include "stretch4.h"
#include "banded.h"

struct Map4 {
    int nt, nd;
    float t0,dt, eps;
    int *x; 
    bool *m;
    float *w, *diag, *offd[3];
    bands slv;
};

map stretch4_init (int n1, float o1, float d1, int nd, float eps)
{
    int i;
    map4 str;
    
    str = (map4) sf_alloc (1, sizeof(*str));

    str->nt = n1; 
    str->t0 = o1; 
    str->dt = d1; 
    str->nd = nd; 
    str->eps = eps;
    
    str->x = sf_intalloc (nd);
    str->m = sf_boolalloc (nd);
    str->w = sf_floatalloc (nd);
    str->diag = sf_floatalloc (n1);
    
    for (i = 0; i < 3; i++) {
	slv->offd[i] = sf_floatalloc (n1-1-i);
    }
  
    str->slv = banded_init (n1,3);

    return str;
}

void stretch4_define (map str, float* coord)
{
    int id, ix, i1, n1;
    float rx, r1;
    const float t=0.01; /* tension */
    
    n1=str->nt;

    for (i1 = 0; i1 < n1-3; i1++) {
	/* regularization */
	str->diag[i1] = str->eps*0.8*(4. - 3.*t);
	str->offd[0][i1] = str->eps*0.15*(11.*t - 12.);
	str->offd[1][i1] = -str->eps*0.24*t;
	str->offd[2][i1] = str->eps*(0.2-0.21*t);
    }
    str->diag[n1-3] = str->eps*0.8*(4. - 3.*t);
    str->offd[0][n1-3] = str->eps*0.15*(11.*t - 12.);
    str->offd[1][n1-3] = -str->eps*0.24*t;
    str->diag[n1-2] = str->eps*0.8*(4. - 3.*t);
    str->offd[0][n1-2] = str->eps*0.15*(11.*t - 12.);
    str->diag[n1-1] = str->eps*0.8*(4. - 3.*t);
    
    for (id = 0; id < str->nd; id++) {
	rx = (coord[id] - str->t0)/str->dt; 
	ix = floorf(rx); 
	rx -= ix;
	if (ix < 1 || ix > str->nt - 2) {
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

void stretch4_apply (map str, float* ord, float* mod)
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

void stretch4_invert (map str, float* ord, float* mod)
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

void stretch4_close (map str)
{
    free (str->x);
    free (str->m);
    free (str->w);
    free (str->diag);
    free (str->offd);
    
    tridiagonal_close (str->slv);
    free (str);
}

/* 	$Id: stretch4.c,v 1.1 2004/04/02 15:55:35 fomels Exp $	 */

