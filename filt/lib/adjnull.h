#ifndef _sf_adjnull_h
#define _sf_adjnull_h

#include "c99.h"

void sf_adjnull (bool adj, bool add, int nx, int ny, float* x, float* y);
void sf_cadjnull (bool adj, bool add, int nx, int ny, 
		  float complex* x, float complex* y);

#endif

/* 	$Id: adjnull.h,v 1.2 2004/03/13 06:00:15 fomels Exp $	 */
