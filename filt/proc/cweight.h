#ifndef _cweight_h
#define _cweight_h

#include <rsf.h>

void cweight_init(float *ww);
void cweight_lop (bool adj, bool add, int nx, int ny, 
		  float complex* xx, float complex* yy);

#endif

/* 	$Id: cweight.h,v 1.1 2004/05/13 22:27:10 fomels Exp $	 */
