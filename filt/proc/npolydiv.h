#ifndef _npolydiv_h
#define _npolydiv_h

#include "nhelix.h"

void npolydiv_init (int nd_in, nfilter aa_in);
void npolydiv_close(void);
void npolydiv_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy);

#endif
