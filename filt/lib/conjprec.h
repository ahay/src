#ifndef _sf_conjprec_h
#define _sf_conjprec_h

#include "c99.h"
#include "_solver.h"

typedef void (*sf_operator2)(int,float*);

void sf_conjprec_init(int nx, int nr, float eps,
		      float tol, bool verb, bool hasp0);
void sf_conjprec_close(void);
void sf_conjprec(sf_operator oper, sf_operator2 prec, 
		 float* p, float* x, const float* dat, int niter);

#endif
