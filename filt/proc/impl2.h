#ifndef _impl2_h
#define _impl2_h

#include <rsf.h>

void impl2_init (float r1, float r2, int n1, int n2, 
		 float tau, float pclip, bool up);
void impl2_close (void);
void impl2_set (float **x);
void impl2_apply (float **x, bool set, bool adj);
void impl2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id: impl2.h,v 1.2 2004/04/09 13:17:10 fomels Exp $	 */
