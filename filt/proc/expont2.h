#ifndef _expont2_h
#define _expont2_h

#include <rsf.h>

void expont2_init(int m1, int m2, float **aa);
void expont2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy);

#endif
