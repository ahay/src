#ifndef _gaussshape2_h
#define _gaussshape2_h

#include <rsf.h>

void gaussshape2_init(int n1, int n2);
void gaussshape2_set(float* a, int n1, int n2, float* pattern, int niter);
void gaussshape2_close(void);
void gaussshape2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id: gaussshape2.h,v 1.1 2004/02/25 16:16:27 fomels Exp $	 */
