#ifndef _interp_h
#define _interp_h

typedef void (*interpolator)(float,int,float*);

void bin_int (float x, int n, float* w);
void lin_int (float x, int n, float* w);
void lg_int (float x, int n, float* w);

#endif

/* 	$Id$	 */
