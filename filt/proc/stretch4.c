#include <math.h>

#include <rsf.h>

#include "stretch4.h"
#include "interp_spline.h"
#include "spline.h"
#include "banded.h"

struct Map4 {
    int nt, nw, nd, ib, ie;
    float t0,dt, eps;
    int *x; 
    bool *m;
    float **w, *diag, *offd[3];
    bands slv;
};

map4 stretch4_init (int n1, float o1, float d1, int nd, int nw, float eps)
{
    int i;
    map4 str;
    
    str = (map4) sf_alloc (1, sizeof(*str));

    str->nt = n1+2*nw; 
    str->t0 = o1; 
    str->dt = d1; 
    str->nd = nd; 
    str->nw = nw;
    str->eps = eps;
    
    str->x = sf_intalloc (nd);
    str->m = sf_boolalloc (nd);
    str->w = sf_floatalloc2 (4,nd);
    str->diag = sf_floatalloc (str->nt);
    
    for (i = 0; i < 3; i++) {
	str->offd[i] = sf_floatalloc (str->nt-1-i);
    }
  
    str->slv = banded_init (str->nt,3);

    return str;
}

void stretch4_define (map4 str, float* coord)
{
    int id, ix, i1, n1, nw, i, j, k;
    float rx, d, o[3], *w;
    
    n1 = str->nt;
    nw = str->nw;

    d = str->eps*2./3.;
    o[0] = -str->eps/8.;
    o[1] = -str->eps/5.;
    o[2] = -str->eps/120.;

    for (i1 = 0; i1 < n1; i1++) {
	/* regularization */
	str->diag[i1] = d;
	for (j=0; j < 3; j++) {
	    if (i1 < n1-1-j) str->offd[j][i1] = o[j];
	}
    }
    
    for (id = 0; id < str->nd; id++) {
	rx = (coord[id] - str->t0)/str->dt; 
	ix = floorf(rx); 
	rx -= ix;
	if (ix < -1 - nw/2 || ix > n1 - nw/2 - 4) {
	    str->m[id] = true; 
	    continue;
	}

	str->x[id] = ix + nw + ; 
	str->m[id] = false; 
	w = str->w[id];

	spline4_int(rx,w);

	k = ix + nw + 1;
	
	for (i=0; i < 4; i++) {
	    str->diag[k+i] += w[i] * w[i];
	    for (j=0; j < 3-i; j++) {
		str->offd[j][k+i] += w[i] * w[i+j+1];
	    }
	}
    }

    banded_define (str->slv, str->diag, str->offd);
    
    str->ib = -2;
    for (i1 = 0; i1 < n1; i1++) {
	if (str->diag[i1] != d) {
	    str->ib = i1-2; 
	    break;
	}
    }
    

    str->ie = n1+2;
    for (i1 = n1-1; i1 >= 0; i1--) {
	if (str->diag[i1] != d) {
	    str->ie = i1+3;
	    break;
	}
    }
}

void stretch4_apply (map4 str, float* ord, float* mod)
{
    int id, it, i, k, nw, nt;
    float *w, *mm;
    
    mm = str->diag;
    nw = str->nw;
    nt = str->nt;

    for (it = 0; it < nt; it++) {
	mm[it] = 0.;
    }
    
    for (id = 0; id < str->nd; id++) {
	if (str->m[id]) continue;
	
	it = str->x[id]; 
	w = str->w[id]; 
	
	for (i=0; i < 4; i++) {
	    k = it + nw + i + 1;
	    mm[k] += w[i]*ord[id];
	}
    }    

    banded_solve (str->slv, mm);

    for (it = 0; it <= str->ib; it++) {
	mm[it] = 0.;
    }
    
    for (it = str->ie; it < nt; it++) {
	mm[it] = 0.;
    }

    spline4_post(nt,nw,nt-nw,mm,mod);
}

void stretch4_close (map4 str)
{
    int i;

    free (str->x);
    free (str->m);
    free (str->w[0]);
    free (str->w);
    free (str->diag);

    for (i = 0; i < 3; i++) {
	free (str->offd[i]);
    }
    
    banded_close (str->slv);
    free (str);
}

/* 	$Id: stretch4.c,v 1.3 2004/04/12 15:40:43 fomels Exp $	 */

