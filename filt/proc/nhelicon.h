#ifndef _nhelicon_h
#define _nhelicon_h

#include "nhelix.h"

void nhelicon_init(nfilter aa);
void nhelicon_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy);

#endif
