#ifndef _cmatmult_h
#define _cmatmult_h

#include <rsf.h>

void cmatmult_init(float complex **bb_in);
void cmatmult_lop (bool adj, bool add, int nx, int ny, 
		   float complex *x, float complex *y);

#endif
