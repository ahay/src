#ifndef _expont_h
#define _expont_h

#include <rsf.h>

void expont_init(int n1_in,int n2_in, float *a1, float *b1);
void expont_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy);

#endif
