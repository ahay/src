#ifndef _mshconest_h
#define _mshconest_h

#include "mshelix.h"
#include <rsf.h>

void mshconest_init(float * x_in, msfilter aa_in);
void mshconest_lop(bool adj, bool add, int na, int ny, float *a, float *y);

#endif

/* 	$Id: mshconest.h,v 1.1 2004/06/11 10:51:33 fomels Exp $	 */
