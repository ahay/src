#ifndef _patching_h
#define _patching_h

#include <rsf.h>

void patching(sf_operator oper, float* modl, float* data, 
	      int dim, int* npatch, int* nwall, int* nwind, float* windwt);

#endif
