#include <rsf.h>

#include "mask.h"

static const bool *m;

void mask_init(const bool *m_in)
{
    m = m_in;
}

void mask_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
{
    int ix;

    if (nx != ny) sf_error("%s: wrong size: %d != %d",nx,ny);

    sf_adjnull (adj,add,nx,ny,x,y);

    for (ix=0; ix < nx; ix++) {
	if (m[ix]) {
	    if (adj) x[ix] += y[ix];
	    else     y[ix] += x[ix];
	}
    }
}

/* 	$Id: mask.c,v 1.4 2004/04/08 14:03:57 fomels Exp $	 */

