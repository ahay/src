#ifndef _int2_h
#define _int2_h

#include <rsf.h>

#include "interp.h"

void  int2_init (float** coord, 
		 float o1, float o2, float d1, float d2,
		 int n1, int n2, 
		 interpolator interp, int nf_in, int nd_in);
void  int2_lop (bool adj, bool add, int nm, int nd, float* x, float* ord);
void int2_close (void);

#endif
