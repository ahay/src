#ifndef _int1_h
#define _int1_h

#include <rsf.h>

#include "interp.h"

void  int1_init (float* coord, float o1, float d1, int n1, 
		 interpolator interp, int nf_in, int nd_in);

void  int1_lop (bool adj, bool add, int nm, int nd, float* x, float* ord);
void int1_close (void);


#endif
