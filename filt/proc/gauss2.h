#ifndef _gauss2_h
#define _gauss2_h

#include <rsf.h>

void gauss2_init(int n1, int n2, float f1, float f2);
void gauss2_close(void);
void gauss2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id: gauss2.h,v 1.1 2004/02/14 06:57:16 fomels Exp $	 */
