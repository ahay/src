#ifndef _gauss_h
#define _gauss_h

#include <rsf.h>

void gauss_init(int n1, float f);
void gauss_close(void);
void gauss_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id: gauss.h,v 1.1 2004/02/14 06:52:41 fomels Exp $	 */
