#include <rsf.h>

#include "matmult.h"

static float** Bb;

/*
  Matmult
  -------
  Simple matrix multiplication operator. 
  Initialized with a pointer to a matrix. */
void matmult_init (float** bb) {
    Bb = bb;
}

void matmult_lop (bool adj, bool add, int nx, int ny, float* x, float*y) {
    int ix, iy;
    sf_adjnull (adj, add, nx, ny, x, y);
    for (ix = 0; ix < nx; ix++) {
	for (iy = 0; iy < ny; iy++) {
	    if (adj) {
		x[ix] += Bb[iy][ix] * y[iy];
	    } else {
		y[iy] += Bb[iy][ix] * x[ix];
	    }
	}
    }
}

