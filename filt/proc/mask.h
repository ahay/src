#ifndef _mask_h
#define _mask_h

#include <rsf.h>

void mask_init(const bool *m_in);
void mask_lop(bool adj, bool add, int nx, int ny, float *x, float *y);

#endif

/* 	$Id: mask.h,v 1.3 2004/04/08 14:03:57 fomels Exp $	 */
