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

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
	if (adj) {
	    xx[i] += yy[i] * w[i];
	} else {
	    yy[i] += xx[i] * w[i];
	}
    }
}

/* 	$Id: weight.c,v 1.2 2004/04/08 14:03:58 fomels Exp $	 */
