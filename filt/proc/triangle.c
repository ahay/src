#include <rsf.h>

#include "triangle.h"

static float *tmp;
static int np, nb, nx;

static void fold (int o, int d, const float *x, float* tmp);
static void doubint (int nx, float *x, bool der);
static void triple (int o, int d, float* x, const float* tmp);

/* 
   Triangle smoothing in 1-D 
   nbox - triangle length
   ndat - data length
*/
void triangle_init (int nbox, int ndat)
{
    nx = ndat;
    nb = nbox;
    np = nx + 2*nb;
    
    tmp = sf_floatalloc(np);
}

static void fold (int o, int d, const float *x, float* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	tmp[i+nb] = x[o+i*d];

    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+(nx-1-i)*d];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+i*d];
    }

    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+i*d];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+(nx-1-i)*d];
    }
}

static void doubint (int nx, float *xx, bool der)
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }

    if (der) return;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }
}

static void triple (int o, int d, float* x, const float* tmp)
{
    int i;
    const float *tmp1, *tmp2;
    float wt;

    tmp1 = tmp + nb;
    tmp2 = tmp + 2*nb;
    wt = 1./(nb*nb);
    
    for (i=0; i < nx; i++) {
	x[o+i*d] = (2. * tmp1[i] - tmp[i] - tmp2[i])*wt;
    }
}

void triangle (int o, int d, bool der, float *x)
{
    fold (o,d,x,tmp);
    doubint (np,tmp,der);
    triple (o,d,x,tmp);
}

void  triangle_close(void)
{
    free (tmp);
}

