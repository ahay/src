#ifndef _loconvol_h
#define _loconvol_h

#include <rsf.h>
#include "helix.h"

void loconvol_init(filter aa);
void loconvol_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy);

#endif
