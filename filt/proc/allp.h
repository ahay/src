#ifndef _allp_h
#define _allp_h

#include <rsf.h>

void allp_init(int nw1, int nj1, int m1, int m2, float **pp1);
void allp_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy);

#endif
