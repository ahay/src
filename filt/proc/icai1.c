#include <rsf.h>

#include "icai1.h"

static int nb, lg;
static float* bb;

/*
  Icai1
  -------
  Internal convolution, adjoint is the input in 1-D. 
  Initialized with the filter and lag, lag=1 is causal. */
void icai1_init (int na, float* aa, int lag) {
    nb = na;
    bb = aa;
    lg = lag;
}

void icai1_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) {
    int b, x, y;

    if(ny != nx) sf_error("%s: size problem: %d != %d",__FILE__,ny,nx);
    sf_adjnull (adj, add, nx, ny, xx, yy);
    
    for( b=0; b < nb; b++) {
	for( y = nb - lg; y <= ny - lg; y++) {
	    x = y - b + lg - 1;
	    if( adj) xx[x] += yy[y] * bb[b];
	    else     yy[y] += xx[x] * bb[b];
	}
    }
}

/* 	$Id: icai1.c,v 1.1 2004/03/26 03:32:17 fomels Exp $	 */
