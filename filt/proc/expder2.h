#ifndef _expder2_h
#define _expder2_h

#include <rsf.h>

void expder2_init(int m1, int m2, float **aa);
void expder2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy);

#endif
