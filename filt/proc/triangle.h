#ifndef _triangle_h
#define _triangle_h

#include <rsf.h>

typedef struct Triangle *triangle;

triangle triangle_init (int nbox, int ndat);
void smooth (triangle tr, int o, int d, bool der, float *x);
void  triangle_close(triangle tr);

#endif

/* 	$Id: triangle.h,v 1.3 2003/10/01 22:45:56 fomels Exp $	 */
