#ifndef _interp_spline_h
#define _interp_spline_h

void spline_int (float x, int n, float* w);
void spline_der (float x, int n, float* w);
void spline4_int (float x, float* w);
void spline4_der (float x, float* w);

#endif

/* 	$Id: interp_spline.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
