#include <rsf.h>

#include "recfilt.h"
#include "adjnull.h"

static float *tt, *aa;
static int na;

/*
  Recfilt
  -------
  Recursive onvolution (polynomial division). 
  Initialized with the filter. */
void recfilt_init( int nd, int nb, float* bb) {
    aa = bb;
    na = nb;
    tt = sf_floatalloc (nd);
}

void recfilt_lop( bool adj, bool add, int nx, int ny, float* xx, float*yy) {
    int ia, iy, ix;

    adjnull( adj, add, nx, ny, xx, yy);

    for (ix=0; ix < nx; ix++) {
	tt[ix] = 0.;
    }

    if (adj) {
	for (ix = nx-1; ix >= 0; ix--) {  
	    tt[ix] = yy[ix];
	    for (ia = 0; ia < na; ia++) {
		iy = ix + ia + 1;
		if( iy >= ny) break;
		tt[ix] -= aa[ia] * tt[iy];
	    }
	}
	for (ix=0; ix < nx; ix++) {
	    xx[ix] += tt[ix];
	}
    } else {
	for (iy = 0; iy < ny; iy++) { 
	    tt[iy] = xx[iy];
	    for (ia = 0; ia < na; ia++) {
		ix = iy - ia - 1;
		if( ix < 0) break;
		tt[iy] -= aa[ia] * tt[ix];
	    }
	}
	for (iy=0; iy < ny; iy++) {
	    yy[iy] += tt[iy];
	}
    }
}

void recfilt_close (void) {
    free (tt);
}

