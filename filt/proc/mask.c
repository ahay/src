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


