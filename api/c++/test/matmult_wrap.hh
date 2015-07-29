#ifndef _matmult_wrap_hh
#define _matmult_wrap_hh

#include <rsf.h>

void matmult_init (int nx, int ny, float** bb);
void matmult_lop (bool adj, bool add, int nx, int ny, float* x, float*y);

#endif
