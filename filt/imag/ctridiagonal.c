#include <math.h>

#include <rsf.h>

#include "ctridiagonal.h"

struct CTris {
    int n;
    float complex *d[2], *o[2], *x[2];  
};

ctris ctridiagonal_init (int n)
{
    ctris slv;
    
    slv = (ctris) sf_alloc (1, sizeof(*slv));
    
    slv->n = n;
    slv->d[0] = sf_complexalloc (n);
    slv->d[1] = sf_complexalloc (n);

    slv->o[0] = sf_complexalloc (n-1);
    slv->o[1] = sf_complexalloc (n-1);

    slv->x[0] = sf_complexalloc (n);
    slv->x[1] = sf_complexalloc (n);

    return slv;
}

void ctridiagonal_define (ctris slv, float complex* diag, float complex* offd)
{
    int k;
    float complex t;

    slv->d[0][0] = diag[0];
    for (k = 1; k < slv->n; k++) {
	t = offd[k-1]; 
	slv->o[0][k-1] = t / slv->d[0][k-1];
	slv->d[0][k] = diag[k] - t * slv->o[0][k-1];
    }
    slv->d[1][slv->n-1] = diag[slv->n-1];
    for (k = slv->n-2; k >= 0; k--) {
	t = offd[k]; 
	slv->o[1][k] = t / slv->d[1][k+1];
	slv->d[1][k] = diag[k] - t * slv->o[1][k];
    }
}

void ctridiagonal_const_define (ctris slv, 
				float complex diag, float complex offd)
{
    int k;
    
    slv->d[0][0] = diag;
    for (k = 1; k < slv->n; k++) {
	slv->o[0][k-1] = offd / slv->d[0][k-1];
	slv->d[0][k] = diag - offd * slv->o[0][k-1];
    }
    slv->d[1][slv->n-1] = diag;
    for (k = slv->n-2; k >= 0; k--) {
	slv->o[1][k] = offd / slv->d[1][k+1];
	slv->d[1][k] = diag - offd * slv->o[1][k];
    }
}

void ctridiagonal_solve (ctris slv, float complex* b)
{
    int k;

    slv->x[0][0] = b[0];
    for (k = 1; k < slv->n; k++) {
	slv->x[0][k] = b[k] - slv->o[0][k-1] * slv->x[0][k-1];
    }
    slv->x[1][slv->n-1] = b[slv->n-1];
    for (k = slv->n-2; k >= 0; k--) {
	slv->x[1][k] = b[k] - slv->o[1][k] * slv->x[1][k+1];
    }
    b[slv->n-1] = slv->x[0][slv->n-1] / slv->d[0][slv->n-1];
    for (k = slv->n-2; k >= slv->n/2; k--) {
	b[k] = slv->x[0][k] / slv->d[0][k] - slv->o[0][k] * b[k+1];
    }
    b[0] = slv->x[1][0] / slv->d[1][0];
    for (k = 1; k < slv->n/2; k++) {
	b[k] = slv->x[1][k] / slv->d[1][k] - slv->o[1][k-1] * b[k-1];
    }
}

void ctridiagonal_close (ctris slv)
{
    free (slv->d[0]); free (slv->d[1]); 
    free (slv->o[0]); free (slv->o[1]); 
    free (slv->x[0]); free (slv->x[1]); 
    free (slv);
}
