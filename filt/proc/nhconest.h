#ifndef _nhconest_h
#define _nhconest_h

#include "nhelix.h"

void nhconest_init(float *x_in, nfilter aa_in, int nhmax_in);
void nhconest_lop(bool adj, bool add, int na, int ny, float *a, float *y);

#endif
