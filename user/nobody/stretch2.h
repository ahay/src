#ifndef _stretch2_h
#define _stretch2_h

#include <rsf.h>

typedef struct Map2 *map2;

map2 stretch2_init (int n1, float o1, float d1, int nd, float eps, float lam);
void stretch2_define (map2 str, float* coord, bool refl);
void stretch2_apply (map2 str, float* ord, float* mod);
void stretch2_invert (map2 str, float* ord, float* mod);
void stretch2_close (map2 str);

#endif

/* 	$Id: stretch2.h,v 1.1 2004/05/25 00:46:12 fomels Exp $	 */
