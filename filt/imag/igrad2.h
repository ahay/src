#ifndef _igrad2_h
#define _igrad2_h

#include <stdbool.h>

void igrad2_init (int n1_in, int n2_in);
void igrad2_lop (bool adj, bool add, int np, int nr, float* p, float* r);

#endif
