#ifndef _neighbors_h
#define _neighbors_h

#include <rsf.h>

void neighbors_init (int *in, float *rdx, int *n, 
		     int order, float *time);
int  neighbours(int i);
int nearsource(float* xs, int* b, float* d, float* vv, bool *plane);

#endif

/* 	$Id: neighbors.h,v 1.3 2004/06/18 01:06:45 fomels Exp $	 */
