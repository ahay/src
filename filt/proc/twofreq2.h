#ifndef _twofreq2_h
#define _twofreq2_h

#include <rsf.h>

void twofreq2_init(int n1, int n2, float eps, float lam, bool gauss);
void twofreq2_close(void);
void twofreq2(int niter, bool verb, float *u, float** pq);

#endif

/* 	$Id: twofreq2.h 704 2004-07-13 18:22:06Z fomels $	 */
