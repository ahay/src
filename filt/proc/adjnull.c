#include <rsf.h>

#include "adjnull.h"

/*
  Function: adjnull
  -----------------
  Zeros out the output (unless add is true). 
  Useful first step for any linear operator. */
void adjnull (bool adj, bool add, int nx, int ny, float* x, float* y) {
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




