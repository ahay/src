#include <rsf.h>

#include "helicon.h"
#include "helix.h"

static filter aa;

/*
  Helicon
  -------
  Helical convolution. 
  Initialized with the filter. */
void helicon_init( filter bb) {
    aa = bb;
}

void helicon_lop( bool adj, bool add, int nx, int ny, float* xx, float*yy) {
    int ia, iy, ix;
    
    sf_adjnull( adj, add, nx, ny, xx, yy);
    
    for (iy = 0; iy < ny; iy++) { /* zero lag */
	if (adj) { 
	    xx[iy] +=	yy[iy];
	} else {
	    yy[iy] +=	xx[iy];
	}
    }

    for (ia = 0; ia < aa->nh; ia++) {
	for (iy = aa->lag[ia]; iy < ny; iy++) {
	    if( aa->mis != NULL && aa->mis[iy]) continue;
	    ix = iy - aa->lag[ia];
	    if(adj) {
		xx[ix] += yy[iy] * aa->flt[ia];
	    } else {
		yy[iy] += xx[ix] * aa->flt[ia];
	    }
	}
    }
}

/* 	$Id: helicon.c,v 1.3 2003/10/21 15:09:08 fomels Exp $	 */
