#ifndef _dip2_h
#define _dip2_h

#include <rsf.h>

void dip2_init(int n1, int n2, float eps, float lam, bool sign, bool gauss);
void dip2_close(void);
void dip2(int niter, int nw, int nj, bool verb, float **u, float** p);

#endif

/* 	$Id: dip2.h,v 1.2 2004/02/26 14:34:25 fomels Exp $	 */
