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

    sf_adjnull (adj,add,nx,ny,x,y);

    for (ix=0; ix < nx; ix++) {
	if (m[ix]) {
	    if (adj) x[ix] += y[ix];
	    else     y[ix] += x[ix];
	}
    }
}

/* 	$Id: mask.c,v 1.3 2003/10/21 15:09:08 fomels Exp $	 */

