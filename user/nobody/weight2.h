#ifndef _weight2_h
#define _weight2_h

#include <rsf.h>

void weight2_init(float *w1, float* w2);
void weight2_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy);

#endif

/* 	$Id: weight2.h,v 1.1 2004/05/25 00:46:12 fomels Exp $	 */
