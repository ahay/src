#include <rsf.h>

#include "helicon.h"
#include "helix.h"
#include "copy.h"

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
    
    copy_lop(adj, add, nx, ny, xx, yy);

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

/* 	$Id: helicon.c,v 1.4 2004/06/18 01:06:45 fomels Exp $	 */
