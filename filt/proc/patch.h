#ifndef _patch_h
#define _patch_h

#include <rsf.h>

void patch_init(int dim_in, int* npatch_in, int* nwall_in, int* nwind_in);
void patch_lop (bool adj, bool add, 
		int nx, int ny, float* wall, float* wind);
void patch_close(void);

#endif
