#include <rsf.h>

#include "polydiv.h"
#include "helix.h"
#include "adjnull.h"

static filter aa;
static float* tt;

/*
  Polydiv
  -------
  Helical inverse convolution (polynomial division). 
  Initialized with the filter. */
void polydiv_init( int nd, filter bb) {
    aa = bb;
    tt = sf_floatalloc (nd);
}

void polydiv_lop( bool adj, bool add, int nx, int ny, float* xx, float*yy) {
    int ia, iy, ix;
    
    adjnull( adj, add, nx, ny, xx, yy);
    
    for (ix=0; ix < nx; ix++) {
	tt[ix] = 0.;
    }

    if (adj) {
	for (ix = nx-1; ix >= 0; ix--) {  
	    tt[ix] = yy[ix];
	    for (ia = 0; ia < aa->nh; ia++) {
		iy = ix + aa->lag[ia];
		if( iy >= ny) continue;
		tt[ix] -= aa->flt[ia] * tt[iy];
	    }
	}
	for (ix=0; ix < nx; ix++) {
	    xx[ix] += tt[ix];
	}
    } else {
	for (iy = 0; iy < ny; iy++) { 
	    tt[iy] = xx[iy];
	    for (ia = 0; ia < aa->nh; ia++) {
		ix = iy - aa->lag[ia];
		if( ix < 0) continue;
		tt[iy] -= aa->flt[ia] * tt[ix];
	    }
	}
	for (iy=0; iy < ny; iy++) {
	    yy[iy] += tt[iy];
	}
    }
}

void polydiv_close (void) {
    free (tt);
}

