#ifndef _sf_conjgrad_h
#define _sf_conjgrad_h

#include "bigsolver.h"
#include "c99.h"

void sf_conjgrad_init(int np, int nx, int nd, int nr, float eps,
		      float tol, bool verb, bool hasp0);
void sf_conjgrad_close(void);

void sf_conjgrad(sf_operator prec, sf_operator oper, sf_operator shape, 
		 float* p, float* x, float* dat, int niter);

#endif

/* 	$Id: conjgrad.h,v 1.3 2004/04/12 15:40:43 fomels Exp $	 */
