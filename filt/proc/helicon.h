#ifndef _helicon_h
#define _helicon_h

#include "helix.h"
#include <rsf.h>

/*
  Helicon
  -------
  Helical convolution. 
  Initialized with the filter. */
void helicon_init( filter aa);
void helicon_lop( bool adj, bool add, int nx, int ny, float* x, float*y);

#endif

/* 	$Id: helicon.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
