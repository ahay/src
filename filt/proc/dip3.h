#ifndef _dip3_h
#define _dip3_h

#include <rsf.h>

void dip3_init(int n1, int n2, int n3, float eps, float lam, bool sign);
void dip3_close(void);
void dip3(int dip, int niter, int nw, int nj, bool verb, 
	  float ***u, float*** p);

#endif

/* 	$Id: dip3.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
