#include <math.h>

#include <rsf.h>

#include "hipass.h"

static float r;

void hipass_init (float eps)
{
    r = 1+0.5*eps*eps;
    r -= sqrtf(r*r-1.);
}

void hipass_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
{
    int i;
    float t;

    sf_adjnull(adj,add,nx,ny,x,y);

    if ( adj) {
	t = y[ny-1];
	x[ny-1] += t;
	for (i=ny-2; i >=0; i--) {
	    t = y[i] - y[i+1] + r*t;
	    x[i] += t;
	}
    } else {
	t = x[0];
	y[0] += t;
	for (i=1; i < nx; i++) {
	    t = x[i] - x[i-1] + r*t;
	    y[i] += t;
	}
    }
}

