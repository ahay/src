#ifndef _trianglen_h
#define _trianglen_h

#include <rsf.h>

void trianglen_init (int dim, int *nbox, int *ndat);
void trianglen_lop (bool adj, bool add, int nx, int ny, float* x, float* y);
void trianglen_close(void);

#endif

/* 	$Id: trianglen.h,v 1.1 2004/04/05 14:38:29 fomels Exp $	 */
