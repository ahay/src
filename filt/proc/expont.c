#include <rsf.h>

#include "expont.h"

static int n1, n2;
static float *a, *b;

void expont_init(int n1_in,int n2_in, float *a1, float *b1)
{
    n1 = n1_in;
    n2 = n2_in;
    a = a1;
    b = b1;
}

void expont_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    int it, ix, i;

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (ix=0; ix < n2; ix++) {
	for (it=2; it < n1; it++) {
	    i = it + ix*n1;
	    if (adj) {
		xx[i-1] -= a[i]*b[i]*yy[i];
		xx[i]   += 0.5*yy[i];
		xx[i-2] += 0.5*b[i]*b[i]*yy[i];
	    } else {
		yy[i] += 0.5*(xx[i] + b[i]*(b[i]*xx[i-2] - 2.*a[i]*xx[i-1]));
	    }
	}
    }
}
