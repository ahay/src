#ifndef _deburst_h
#define _deburst_h

#include <rsf.h>

void deburst (int n, int niter, sf_weight wght, float eps, 
	      const float *dd, float *hh);

#endif
