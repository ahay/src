#ifndef _boxfilter_h
#define _boxfilter_h

#include "helix.h"

void box (int dim, const int *nd, const int *center, const int *na, 
	  const filter aa, int nc, float* cube);

#endif
