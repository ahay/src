#ifndef _dip2_h
#define _dip2_h

#include <rsf.h>

void dip2_init(int n1, int n2, float eps, float lam, bool sign, bool gauss);
void dip2_close(void);
void dip2(int niter, int nw, int nj, bool verb, 
	  float **u, float** p, bool **mask);

#endif

/* 	$Id$	 */
