#ifndef _npolydiv2_h
#define _npolydiv2_h

#include "nhelix.h"

void npolydiv2_init (int nd_in, nfilter aa_in);
void npolydiv2_close(void);
void npolydiv2_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy);

#endif
