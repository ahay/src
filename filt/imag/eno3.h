#ifndef _eno3_h
#define _eno3_h

#include "eno.h"

typedef struct Eno3 *eno3;

eno3 eno3_init (int order, int n1, int n2, int n3);
void eno3_set (eno3 pnt, float*** c);
void eno3_set1 (eno3 pnt, float* c);
void eno3_close (eno3 pnt);
void eno3_apply (eno3 pnt, int i, int j, int k,
		 float x, float y, float z,
		 float* f, float* f1, der what);

#endif

/* 	$Id: eno3.h,v 1.2 2003/09/30 14:30:52 fomels Exp $	 */
