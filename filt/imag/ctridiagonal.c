#include <math.h>

#include <rsf.h>

#include "ctridiagonal.h"

struct CTris {
    int n;
    float complex *d, *o;  
};

ctris ctridiagonal_init (int n)
{
    ctris slv;
    
    slv = (ctris) sf_alloc (1, sizeof(*slv));
    
    slv->n = n;
    slv->d = sf_complexalloc (n);
    slv->o = sf_complexalloc (n-1);

    return slv;
}

void ctridiagonal_define (ctris slv, float complex* diag, float complex* offd)
{
    int k;
    float complex t;

    slv->d[0] = diag[0];
    for (k = 1; k < slv->n; k++) {
	t = offd[k-1]; 
	slv->o[k-1] = t / slv->d[k-1];
	slv->d[k] = diag[k] - t * slv->o[k-1];
    }
}

void ctridiagonal_const_define (ctris slv, 
				float complex diag, float complex offd)
{
    int k;
    
    slv->d[0] = diag;
    for (k = 1; k < slv->n; k++) {
	slv->o[k-1] = offd / slv->d[k-1];
	slv->d[k] = diag - offd * slv->o[k-1];
    }
}

void ctridiagonal_solve (ctris slv, float complex* b)
{
    int k, n;

    n = slv->n;

    for (k = 1; k < n; k++) {
	b[k] -= slv->o[k-1] * b[k-1];
    }
    b[n-1] /= slv->d[n-1];
    for (k = n-2; k >= 0; k--) {
	b[k] = b[k] / slv->d[k] - slv->o[k] * b[k+1];
    }
}

void ctridiagonal_close (ctris slv)
{
    free (slv->d); 
    free (slv->o); 
    free (slv);
}

/* 	$Id: ctridiagonal.c,v 1.4 2003/09/30 14:30:52 fomels Exp $	 */
