#ifndef _peftc_h
#define _peftc_h

#include <rsf.h>

void peftc_init (int na, int ny, float *aa, float *yy);
void peftc_lop (bool adj, bool add, int nx, int nr, float *x, float *r);

#endif
