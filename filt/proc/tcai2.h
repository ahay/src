#ifndef _tcai2_h
#define _tcai2_h

#include <rsf.h>

void tcai2_init (int na, int nx, float *aa);
void tcai2_lop (bool adj, bool add, int nx, int nr, float *x, float *r);

#endif
