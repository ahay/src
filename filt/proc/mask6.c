#include <rsf.h>

#include "mask6.h"

void mask6_apply (int nw, int nj1, int nj2, int nx, int ny, float **yy) 
{
    int ix, iy, iw, is;
    float **xx;

    xx = sf_floatalloc2(nx,ny);
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    if (yy[iy][ix] == 0.) {
		yy[iy][ix] = 1.;
	    } else {
		yy[iy][ix] = 0.;
	    }
	    xx[iy][ix] = 0.;
	}
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj1;
		xx[iy][ix] += (yy[iy+1][ix+is] + 
			       yy[iy  ][ix-is]);
	    }
	}
    }
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    yy[iy][ix] = 0.;
	}
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj2;
		yy[iy][ix] += (xx[iy+1][ix+is] + 
			       xx[iy  ][ix-is]);
	    }
	}
    }

    free(xx[0]);
    free(xx);
}
