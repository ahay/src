#ifndef _sf_cconjgrad_h
#define _sf_cconjgrad_h

#include "bigsolver.h"
#include "c99.h"

void sf_cconjgrad_init(int np, int nx, int nd, int nr, float eps,
		      float tol, bool verb, bool hasp0);
void sf_cconjgrad_close(void);

void sf_cconjgrad(sf_coperator prec, sf_coperator oper, sf_coperator shape, 
		  float complex* p, float complex* x, float complex* dat, 
		  int niter);

#endif

/* 	$Id: cconjgrad.h,v 1.1 2004/05/13 22:26:56 fomels Exp $	 */
