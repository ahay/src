#ifndef _spline_h
#define _spline_h

#include "banded.h"
#include "tridiagonal.h"

bands spline_init (int nw, int nd);
tris spline4_init (int nd);
void spline4_post (int n, int n1, int n2, const float* inp, float* out);

void spline2 (bands slv1, bands slv2, int n1, int n2, float** dat, float* tmp);

#endif

/* 	$Id: spline.h,v 1.3 2004/04/03 02:41:17 fomels Exp $	 */
