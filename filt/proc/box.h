#ifndef _box_h
#define _box_h

#include <rsf.h>

void box_init (int nbox, int ndat, bool lin);
void boxsmooth  (int o, int d, float *x, float *y);
void boxsmooth2 (int o, int d, float *x, float *y);
void box_lop(bool adj, bool add, int nx, int ny, float* x, float* y);
void box_close (void);

#endif

/* 	$Id: box.h,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */
