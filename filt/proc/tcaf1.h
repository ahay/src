#ifndef _tcaf1_h
#define _tcaf1_h

#include <rsf.h>

void tcaf1_init(int ny, float* yy);
void tcaf1_lop(bool adj, bool add, int nb, int ny, float *bb, float *yy);

#endif
