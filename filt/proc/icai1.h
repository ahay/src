#ifndef _icai1_h
#define _icai1_h

#include <rsf.h>

/*
  Icai1
  -------
  Internal convolution, adjoint is the input in 1-D. 
  Initialized with the filter and the lag, lag=1 is causal. */
void icai1_init( int na, float* aa, int lag);
void icai1_lop( bool adj, bool add, int nx, int ny, float* x, float*y);

#endif

/* 	$Id: icai1.h,v 1.1 2004/03/26 03:32:17 fomels Exp $	 */
