#ifndef _tcai1_h
#define _tcai1_h

#include <rsf.h>

/*
  Tcai1
  -------
  Transient convolution, adjoint is the input in 1-D. 
  Initialized with the filter. */
void tcai1_init( int na, float* aa);
void tcai1_lop( bool adj, bool add, int nx, int ny, float* x, float*y);

#endif
