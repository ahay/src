#ifndef _hconest_h
#define _hconest_h

#include "helix.h"
#include <rsf.h>

void hconest_init(float * x_in, filter aa_in);
void hconest_lop(bool adj, bool add, int na, int ny, float *a, float *y);

#endif
