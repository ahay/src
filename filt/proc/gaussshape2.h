#ifndef _gaussshape2_h
#define _gaussshape2_h

#include <rsf.h>

void gaussshape2_init(int n1, int n2);
void gaussshape2_set(float* a, const float* pattern, int niter);
void gaussshape2_set2(const float* a);
void gaussshape2_close(void);

#endif

/* 	$Id: gaussshape2.h,v 1.2 2004/02/26 05:16:08 fomels Exp $	 */
