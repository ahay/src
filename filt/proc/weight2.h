#ifndef _weight2_h
#define _weight2_h

#include <rsf.h>

void weight2_init(int nw1, int n, float *ww);
void weight2_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy);
void weight2_close(void);

#endif

/* 	$Id$	 */
