#ifndef _polydiv_h
#define _polydiv_h

#include "helix.h"
#include <rsf.h>

/*
  Polydiv
  -------
  Helical inverse convolution (polynomial division). 
  Initialized with the filter. */
void polydiv_init( int nd, filter aa);
void polydiv_lop( bool adj, bool add, int nx, int ny, float* x, float*y);
void polydiv_close( void);

#endif

/* 	$Id: polydiv.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
