#include <assert.h>

#include <rsf.h>

#include "weight.h"

static float* w;

void weight_init(float *w1)
{
    w = w1;
}

void weight_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
{
    int i;

    assert(ny==nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
	if (adj) {
	    xx[i] += yy[i] * w[i];
	} else {
	    yy[i] += xx[i] * w[i];
	}
    }
}

/* 	$Id: weight.c,v 1.1 2004/01/15 02:38:44 fomels Exp $	 */
