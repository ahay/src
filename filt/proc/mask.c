#include <rsf.h>

#include "mask.h"
#include "adjnull.h"

static const bool *m;

void mask_init(const bool *m_in)
{
    m = m_in;
}

void mask_lop(bool adj, bool add, int nx, int ny, float *x, float *y)
{
    int ix;

    adjnull (adj,add,nx,ny,x,y);

    for (ix=0; ix < nx; ix++) {
	if (m[ix]) {
	    if (adj) x[ix] += y[ix];
	    else     y[ix] += x[ix];
	}
    }
}

/* 	$Id: mask.c,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */

