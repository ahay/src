#ifndef _neighbors_h
#define _neighbors_h

void neighbors_init (int *in, float *rdx, int *n, 
		     int order, float *time);
int  neighbours(int i);
int nearsource(float* xs, int* b, float* d, float* vv);

#endif

/* 	$Id: neighbors.h,v 1.2 2003/09/30 14:30:53 fomels Exp $	 */
