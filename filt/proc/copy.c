#include <rsf.h>

#include "copy.h"

void copy_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
{
    int i;

    if (ny!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,ny,nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    for (i=0; i < nx; i++) {
	if (adj) {
	    xx[i] += yy[i];
	} else {
	    yy[i] += xx[i];
	}
    }
}

/* 	$Id: copy.c,v 1.1 2004/05/18 11:41:27 fomels Exp $	 */
