#include "adjnull.h"

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

/* 	$Id: adjnull.c,v 1.2 2004/03/13 06:00:15 fomels Exp $	 */

