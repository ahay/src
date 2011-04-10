#include <rsf.h>

#include "scale.h"

static float w;

void scale_init(float w1)
{
    w = w1;
}

void scale_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
	if (adj) {
	    xx[i] += yy[i] * w;
	} else {
	    yy[i] += xx[i] * w;
	}
    }
}

/* 	$Id$	 */
