#include <rsf.h>

#include "cmatmult.h"  

static float complex **bb;

void cmatmult_init(float complex **bb_in)
{
    bb = bb_in;
}

void cmatmult_lop (bool adj, bool add, int nx, int ny, 
		   float complex *x, float complex *y)
{
    int ix, iy;
    sf_cadjnull (adj,add,nx,ny,x,y);

    for (ix=0; ix < nx; ix++) {
	for (iy=0; iy < ny; iy++) {
	    if (adj) {
		x[ix] +=conjf(bb[iy][ix])*y[iy];
	    } else {
		y[iy] += bb[iy][ix]*x[ix];
	    }
	}
    }
}

