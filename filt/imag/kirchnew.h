#ifndef _kirchnew_h
#define _kirchnew_h

#include <rsf.h>

void kirchnew_init (float *vrms_in, float t0_in, float dt_in, float dx_in, 
		    int nt_in, int nx_in, int sw_in);
void kirchnew_lop (bool adj, bool add, int nm, int nd, 
		   float *modl, float *data);

#endif
