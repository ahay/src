#ifndef _trisl_h
#define _trisl_h

#include <rsf.h>

void trisl_init(int m1, int m2, int rect1, int rect2);
void trisl_set (float** p);
void trisl_close(void);
void trisl_lop(bool adj, bool add, int nx, int ny, float* x, float* y);

#endif
