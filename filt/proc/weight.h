#ifndef _weight_h
#define _weight_h

#include <rsf.h>

void weight_init(float *ww);
void weight_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy);

#endif

/* 	$Id: weight.h,v 1.2 2004/02/27 20:59:57 fomels Exp $	 */
