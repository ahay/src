#include <assert.h>

#include <rsf.h>

#include "weight2.h"

static float *w1, *w2;

void weight2_init(float *ww1, float *ww2)
{
    w1 = ww1;
    w2 = ww2;
}

void weight2_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
{
    int i;

    if (2*ny != nx) sf_error("%s: size mismatch: 2*%d != %d",__FILE__,ny,nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < ny; i++) {
	if (adj) {
	    xx[i]    += yy[i] * w1[i];
	    xx[i+ny] += yy[i] * w2[i];
	} else {
	    yy[i] += xx[i] * w1[i] + xx[i+ny] * w2[i];
	}
    }
}

/* 	$Id: weight2.c,v 1.1 2004/05/25 00:46:12 fomels Exp $	 */
