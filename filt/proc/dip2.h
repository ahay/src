#ifndef _dip2_h
#define _dip2_h

#include <rsf.h>

void dip2_init(int n1, int n2, float eps, float lam, bool sign);
void dip2_close(void);
void dip2(int niter, int nw, int nj, bool verb, float **u, float** p);

#endif

/* 	$Id: dip2.h,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */
