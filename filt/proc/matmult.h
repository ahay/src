#ifndef _matmult_h
#define _matmult_h

#include <rsf.h>

/*
  Matmult
  -------
  Simple matrix multiplication operator. 
  Initialized with a pointer to a matrix. */
void matmult_init (float** bb);
void matmult_lop (bool adj, bool add, int nx, int ny, float* x, float*y);

#endif
