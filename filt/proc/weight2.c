#include <assert.h>

#include <rsf.h>

#include "weight2.h"

static int nw;
static float **w;

void weight2_init(int nw1, int n, float *ww)
{
    int iw;

    nw = nw1;
    w = (float**) sf_alloc(nw,sizeof(float*));

    for (iw=0; iw < nw; iw++) {
	w[iw] = ww+iw*n;
    }
}

void weight2_close(void)
{
    free(w);
}

void weight2_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
{
    int i, iw;

    if (nw*ny != nx) sf_error("%s: size mismatch: %d*%d != %d",
			      __FILE__,nw,ny,nx);

    sf_adjnull (adj, add, nx, ny, xx, yy);
  
    for (iw=0; iw < nw; iw++) {
	for (i=0; i < ny; i++) {
	    if (adj) {
		xx[i+iw*ny] += yy[i] * w[iw][i];
	    } else {
		yy[i] += xx[i+iw*ny] * w[iw][i];
	    }
	}
    }
}

/* 	$Id$	 */
