#include <rsf.h>

#include "box.h"

static  int np, nb, nx;
static float *tmp;

static void causint (int nx, float *x);
static void causint2 (int nx, float *x);
static void duble (int o, int d, int nx, int nb, 
		    float* x, const float* tmp);
static void duble2 (int o, int d, int nx, int nb, const float* x, float* tmp);

void box_init (int nbox, int ndat, bool lin)
{
    nx = ndat;
    nb = nbox;
    np = ndat + nbox;
    
    if (lin) tmp = sf_floatalloc(np);
}

void box_close (void)
{
    free (tmp);
}

static void causint (int nx, float *xx)
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }
}

static void causint2 (int nx, float *xx)
{
    int i;
    float t;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }
}

static void duble (int o, int d, int nx, int nb, float* x, const float* tmp)
{
    int i;
    const float *tmp1;
    float wt;

    tmp1 = tmp + nb;
    wt = 1./nb;
    
    for (i=0; i < nx; i++) {
	x[o+i*d] = (tmp[i] - tmp1[i])*wt;
    }
}

static void duble2 (int o, int d, int nx, int nb, const float* x, float* tmp)
{
    int i;
    float *tmp1;
    float wt;

    tmp1 = tmp + nb;
    wt = 1./nb;
    
    for (i=0; i < nx + nb; i++) {
	tmp[i] = 0;
    }

    for (i=0; i < nx; i++) {
	tmp1[i] -= x[o+i*d]*wt;
	tmp[i]  += x[o+i*d]*wt; 
    }
}

void boxsmooth (int o, int d, float *x, float *y)
{
    causint (np,y);
    duble (o,d,nx,nb,x,y);
}

void boxsmooth2 (int o, int d, float *x, float *y)
{
    duble2 (o,d,nx,nb,x,y);
    causint2 (np,y);
}

void box_lop(bool adj, bool add, int nx, int ny, float* x, float* y) 
{
    int i;

    sf_adjnull(adj,add,nx,ny,x,y);

    if (adj) {
	boxsmooth(0,1,tmp,y);
	for (i=0; i < nx; i++) {
	    x[i] += tmp[i];
	}
    } else {
	boxsmooth2(0,1,x,tmp);
	for (i=0; i < ny; i++) {
	    y[i] += tmp[i];
	}
    }
}

/* 	$Id: box.c,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */

