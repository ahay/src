#include <rsf.h>

#include "tcai1.h"
#include "adjnull.h"

static int nb;
static float* bb;

/*
  Tcai1
  -------
  Transient convolution, adjoint is the input in 1-D. 
  Initialized with the filter. */
void tcai1_init (int na, float* aa) {
    nb = na;
    bb = aa;
}

void tcai1_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) {
    int b, x, y;

    if(ny < nx+nb) sf_error("%s: size problem: %d < %d+%d",
			    __FILE__,ny,nx,nb);
    adjnull (adj, add, nx, ny, xx, yy);
    
    for( b=0; b < nb; b++) {
	for( x=0; x < nx; x++) {
	    y = x + b;
	    if( adj) xx[x] += yy[y] * bb[b];
	    else     yy[y] += xx[x] * bb[b];
	}
    }
}

