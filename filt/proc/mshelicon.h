#ifndef _mshelicon_h
#define _mshelicon_h

#include "mshelix.h"
#include <rsf.h>

/*
  mshelicon
  -------
  multi-scale helical convolution. 
  Initialized with the filter. */
void mshelicon_init( msfilter aa);
void mshelicon_lop( bool adj, bool add, int nx, int ny, float* x, float*y);

#endif

/* 	$Id: mshelicon.h,v 1.1 2004/06/11 10:51:34 fomels Exp $	 */
