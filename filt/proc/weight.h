#ifndef _weight_h
#define _weight_h

#include <rsf.h>

void weight_init(float *ww);
void weight_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy);

#endif
