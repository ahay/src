#include <math.h>

#include <rsf.h>

#include "steepdip.h"
#include "helix.h"

filter steep(int dim, int *n, int *a, float *d, float vel, float tgap) 
{
    int *lag, c[SF_MAX_DIM], i, h, na, nx, it, j;
    float x, t0;
    filter aa;

    na = 1;
    for (j=0; j < dim; j++) {
	na *= a[j];
    }
    nx = dim-1;

    it = dim; 
    lag = sf_intalloc(na);

    for (h=i=0; i < na; i++) { 
	sf_line2cart(dim, a, i, c);

	for (j=0; j < dim-1; j++) {
	    c[j] -= (a[j]+1)/2;
	}

	t0 = 0.;
	for (j=0; j < dim; j++) {
	    x = d[j]*c[j];
	    t0 += x*x;
	}
	t0 = sqrtf(t0)/vel;
	if (t0 < tgap) t0 = tgap;

	if(x > t0) continue;

	lag[h++] = sf_cart2line(dim, n, c);
    }

    aa = allocatehelix(h);

    for (i=0; i < h; i++) {
	aa->lag[i] = lag[i]-1;
	aa->flt[i] = 0.;
    }

    free(lag);

    return aa;
}

