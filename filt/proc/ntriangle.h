#ifndef _triangle_h
#define _triangle_h

#include <rsf.h>

typedef struct NTriangle *ntriangle;

ntriangle ntriangle_init (int nbox, int ndat);
void nsmooth  (ntriangle tr, int o, int d, bool der, const int *t, float *x);
void nsmooth2 (ntriangle tr, int o, int d, bool der, const int *t, float *x);
void  ntriangle_close(ntriangle tr);

#endif

/* 	$Id: ntriangle.h 691 2004-07-04 19:28:08Z fomels $	 */
