#ifndef _regrid_h
#define _regrid_h

#include "helix.h"

/* Regrid
   ------
   Convert a helix filter from one data size to another. */
void regrid( int dim, const int* nold, const int* nnew, filter aa);

#endif

/* 	$Id: regrid.h,v 1.3 2003/10/01 22:45:56 fomels Exp $	 */
