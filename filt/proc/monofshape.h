#ifndef _monofshape_h
#define _monofshape_h

#include <rsf.h>

void monofshape_init(int n1);
void monofshape_set (float f, int n, float* pattern, int niter);
void monofshape_close(void);
void monofshape_lop (bool adj, bool add, int nx, int ny, float* x, float* y);

#endif

/* 	$Id: monofshape.h,v 1.1 2004/02/14 06:57:16 fomels Exp $	 */
