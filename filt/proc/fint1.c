#include <rsf.h>

#include "fint1.h"
#include "extend.h"
#include "spline.h"
#include "interp_spline.h"
#include "tridiagonal.h"

struct Vint1 {
    float** spl;
    float w[4];
    int n1, nw, dim;
    tris slv;
};

struct Fint1 {
    float *spl, w[4];
    int n1, nw;
    tris slv;
};

float* fint1_coeff (fint1 fnt, int n)
{
    float* coeff;
    coeff = &(fnt->spl[n+fnt->nw]);

    return coeff;
}

float* vint1_coeff (vint1 fnt, int n, int dim)
{
    float* coeff;
    coeff = &(fnt->spl[dim][n+fnt->nw]);
    
    return coeff;
}

vint1 vint1_init (int nw, int n1, int dim)
{
    vint1 fnt;

    fnt = (vint1) sf_alloc (1, sizeof (*fnt));

    fnt->nw = nw; 
    fnt->n1 = n1; 
    fnt->dim = dim;
    fnt->spl = sf_floatalloc2 (n1+2*nw,dim);
    fnt->slv = spline4_init (n1+2*nw);
    
    return fnt;
}

fint1 fint1_init (int nw, int n1)
{
    fint1 fnt;
    
    fnt = (fint1) sf_alloc (1, sizeof (*fnt));
    
    fnt->nw = nw; 
    fnt->n1 = n1; 
    fnt->spl = sf_floatalloc (n1+2*nw);
    fnt->slv = spline4_init (n1+2*nw);
    
    return fnt;
}

void vint1_set (vint1 fnt, float** dat)
{
    int i;
    for (i = 0; i < fnt->dim; i++) {
	extend (fnt->nw,fnt->n1,dat[i],fnt->spl[i]);
	fnt->spl[i][0] *= (5./6.);
	fnt->spl[i][fnt->n1+2*fnt->nw-1] *= (5./6.);
	tridiagonal_solve (fnt->slv,fnt->spl[i]);
    }
}

void fint1_set (fint1 fnt, float* dat)
{
    extend (fnt->nw,fnt->n1,dat,fnt->spl);
    fnt->spl[0] *= (5./6.);
    fnt->spl[fnt->n1+2*fnt->nw-1] *= (5./6.);
    tridiagonal_solve (fnt->slv,fnt->spl);
}

void fint1_close (fint1 fnt)
{
    free (fnt->spl);
    tridiagonal_close (fnt->slv);
    free (fnt);
}

void vint1_close (vint1 fnt)
{
    free (fnt->spl[0]);
    free (fnt->spl);
    tridiagonal_close (fnt->slv);
    free (fnt);
}

float fint1_apply (fint1 fnt, int i, float x, bool der) 
{
    float f;
    int j, k;

    if (der) {
	spline4_der (x,fnt->w);
    } else {
	spline4_int (x,fnt->w);
    }
    
    f = 0.;
    for (j = 0; j < 4; j++) {
	k = i+fnt->nw/2+j+1;
	f += fnt->w[j]*fnt->spl[k];
    }
    return f;
}

void vint1_apply (vint1 fnt, int i, float x, bool der, float* f) 
{
    int j, k, n;
    
    if (der) {
	spline4_der (x,fnt->w);
    } else {
	spline4_int (x,fnt->w);
    }
    
    for (n=0; n < fnt->dim; n++) {
	f[n] = 0.;
	for (j = 0; j < 4; j++) {
	    k = i+fnt->nw/2+j+1;
	    f[n] += fnt->w[j]*fnt->spl[n][k];
	}
    }
}

void stretch(fint1 str, 
	     float (*map)(float),
	     int n1, float d1, float o1,
	     int n2, float d2, float o2,
	     float *trace)
{
    int i2, it;
    float t;

    for (i2=0; i2 < n2; i2++) {
	t = o2+i2*d2;
	t = map(t);
	t = (t-o1)/d1;
	it = t;
	if (it >= 0 && it < n1) {
	    trace[i2] = fint1_apply(str,it,t-it,false);
	} else {
	    trace[i2] = 0.;
	}
    }
}

/* 	$Id: fint1.c,v 1.3 2004/04/01 02:12:42 fomels Exp $	 */
