#include <rsf.h>

#include "cweight.h"

static float* w;

void cweight_init(float *w1)
{
    w = w1;
}

void cweight_lop (bool adj, bool add, int nx, int ny, 
		 float complex* xx, float complex* yy)
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
	if (adj) {
	    xx[i] += yy[i] * w[i];
	} else {
	    yy[i] += xx[i] * w[i];
	}
    }
}

/* 	$Id: cweight.c,v 1.1 2004/05/13 22:27:10 fomels Exp $	 */
