#ifndef _spline_h
#define _spline_h

#include "banded.h"
#include "tridiagonal.h"

bands spline_init (int nw, int nd);
tris spline4_init (int nd);

void spline2 (bands slv1, bands slv2, int n1, int n2, float** dat, float* tmp);

#endif

/* 	$Id: spline.h,v 1.2 2003/10/01 22:45:56 fomels Exp $	 */
