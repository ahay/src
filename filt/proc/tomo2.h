#ifndef _tomo2_h
#define _tomo2_h

#include "int2.h"

void tomo2_init (float*** rays, int *raylen, int nrays, 
		 float o1, float o2, float d1, float d2,
		 int n1, int n2, 
		 interpolator interp, int nf_in);
void  tomo2_lop (bool adj, bool add, int nm, int ny, float* x, float* ord);
void tomo2_close (void);

#endif

/* 	$Id: tomo2.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
