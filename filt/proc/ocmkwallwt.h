#ifndef _ocmkwallwt_h
#define _ocmkwallwt_h

#include <stdio.h>

#include <rsf.h>

void ocmkwallwt(bool inv, int dim, int* npatch, int* nwall, int* nwind, 
		float* windwt, FILE* wallwt);

#endif
