#ifndef _scale_h
#define _scale_h

#include <rsf.h>

void scale_init(float ww);
void scale_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy);

#endif

/* 	$Id: scale.h,v 1.1 2004/05/25 00:46:12 fomels Exp $	 */
