#ifndef _dip3_h
#define _dip3_h

#include <rsf.h>

void dip3_init(int n1, int n2, int n3, int* rect, int niter, bool sign1);
void dip3_close(void);
void dip3(int dip, int niter, int nw, int nj, bool verb, 
	  float ***u, float*** p, bool*** m);

#endif

/* 	$Id$	 */
