#include "adjnull.h"

#include "c99.h"

/*
  Function: adjnull
  -----------------
  Zeros out the output (unless add is true). 
  Useful first step for any linear operator. */
void sf_adjnull (bool adj, bool add, int nx, int ny, float* x, float* y) {
    int i;
    
    if(add) return;
    
    if(adj) {
	for (i = 0; i < nx; i++) {
	    x[i] = 0.;
	}
    } else {
	for (i = 0; i < ny; i++) {
	    y[i] = 0.;
	}
    }
}

/* test */
void sf_cadjnull (bool adj, bool add, int nx, int ny, 
		  float complex* x, float complex* y) {
    int i;
    
    if(add) return;
    
    if(adj) {
	for (i = 0; i < nx; i++) {
	    x[i] = 0.;
	}
    } else {
	for (i = 0; i < ny; i++) {
	    y[i] = 0.;
	}
    }
}

/* 	$Id: adjnull.c,v 1.3 2004/04/02 02:22:38 fomels Exp $	 */

