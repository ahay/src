#ifndef _monof2_h
#define _monof2_h

#include <rsf.h>

void monof2(float **data, int niter, float* a, 
	    int nx, float dx, float x0, 
	    int ny, float dy, float y0, 
	    bool verb);

#endif

/* 	$Id: monof2.h,v 1.1 2004/02/25 16:16:27 fomels Exp $	 */
