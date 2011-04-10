#include <math.h>

#include <rsf.h>

#include "stretch2.h"
#include "banded.h"

struct Map2 {
    int nt, nd;
    float t0,dt, eps, lam;
    int *x; 
    bool *m;
    float *w, *diag, **offd;
    bands slv;
};

map2 stretch2_init (int n1, float o1, float d1, int nd, float eps, float lam)
{
    map2 str;
    
    str = (map2) sf_alloc (1, sizeof(*str));

    str->nt = n1; 
    str->t0 = o1; 
    str->dt = d1; 
    str->nd = nd; 
    str->eps = eps;
    str->lam = lam;
    
    str->x = sf_intalloc (nd);
    str->m = sf_boolalloc (nd);
    str->w = sf_floatalloc (nd);
    str->diag = sf_floatalloc (n1);
    str->offd = sf_floatalloc2 (n1,2);
  
    str->slv = banded_init (n1,2);

    return str;
}

void stretch2_define (map2 str, float* coord, bool refl)
{
    int id, ix, i1;
    float rx, r1;
    
    for (i1 = 0; i1 < str->nt; i1++) {
	/* regularization */
	str->diag[i1] = str->eps*(6.*str->lam+2.5*(1.-str->lam));
	str->offd[0][i1] = -4.*str->eps*(str->lam+(1.-str->lam)/3.);
	str->offd[1][i1] = str->eps*(str->lam+(1.-str->lam)/12.);
    }

    if (refl) {
	sf_warning("reflecting boundary conditions");
	str->diag[0] = str->diag[str->nt-1] = 
	    str->eps*(2.*str->lam+7.*(1.-str->lam)/6.);
	str->offd[0][0] = str->offd[0][str->nt-2] = 
	    -str->eps*(3.*str->lam+1.25*(1.-str->lam));
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
	
	str->diag[ix]    += r1 * r1;
	str->diag[ix+1]  += rx * rx;
	str->offd[0][ix] += r1 * rx;
    }
    
    banded_define (str->slv, str->diag, str->offd);
}

void stretch2_apply (map2 str, float* ord, float* mod)
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
    
    banded_solve (str->slv, mod);
}

void stretch2_invert (map2 str, float* ord, float* mod)
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

void stretch2_close (map2 str)
{
    free (str->x);
    free (str->m);
    free (str->w);
    free (str->diag);
    free (str->offd[0]);
    free (str->offd);
    
    banded_close (str->slv);
    free (str);
}

/* 	$Id$	 */

