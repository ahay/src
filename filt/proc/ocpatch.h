#ifndef _ocpatch_h
#define _ocpatch_h

#include <stdio.h>

#include <rsf.h>

void ocpatch_init(int dim, int nw, int np, 
		  int* npatch, int* nwall, int* nwind);
void ocpatch_close(void);
void ocpatch_lop (int ip, bool adj, FILE *wall, float* wind);
void ocpatch_flop (int ip, bool adj, sf_file wall, float* wind);

#endif
