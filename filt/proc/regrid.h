#ifndef _regrid_h
#define _regrid_h

#include "helix.h"

/* Regrid
   ------
   Convert a helix filter from one data size to another. */
void regrid( int dim, int* nold, int* nnew, filter aa);

#endif
