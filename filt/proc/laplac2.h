#ifndef _laplac2_h
#define _laplac2_h

#include <rsf.h>

void laplac2_init(int m1, int m2);
void laplac2_lop(bool adj, bool add, int np, int nr, float *p, float *r);

#endif
