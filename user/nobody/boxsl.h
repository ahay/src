#ifndef _boxsl_h
#define _boxsl_h

#include <rsf.h>

void boxsl_init(int m1, int m2, int rect1, int rect2);
void boxsl_set (int m2, float** p);
void boxsl_close(void);
void boxsl_lop(bool adj, bool add, int nx, int ny, float* x, float* y);
void trisl_lop(bool adj, bool add, int nx, int ny, float* x, float* y);

#endif
