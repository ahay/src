#ifndef _twodip2_h
#define _twodip2_h

#include <rsf.h>

void twodip2_init(int n1, int n2, float eps, float lam, bool sign, bool gauss);
void twodip2_close(void);
void twodip2(int niter, int nw, int nj1, int nj2, 
	     bool verb, float **u, float*** pq);

#endif

/* 	$Id: twodip2.h,v 1.1 2004/05/25 00:46:12 fomels Exp $	 */
