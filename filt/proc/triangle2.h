#ifndef _triangle2_h
#define _triangle2_h

#include <rsf.h>

void triangle2_init (int nbox1, int nbox2, int ndat1, int ndat2);
void twotriangle2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);
void triangle2_lop (bool adj, bool add, int nx, int ny, float* x, float* y);
void triangle2_close(void);

#endif

/* 	$Id: triangle2.h,v 1.3 2004/02/27 20:59:57 fomels Exp $	 */
