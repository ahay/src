#ifndef _sf_adjnull_h
#define _sf_adjnull_h

#include "c99.h"

void sf_adjnull (bool adj, bool add, int nx, int ny, float* x, float* y);

#ifndef __cplusplus

void sf_cadjnull (bool adj, bool add, int nx, int ny, 
		  float complex* x, float complex* y);

#endif /* c++ */

#endif

/* 	$Id: adjnull.h,v 1.3 2004/06/23 08:54:31 fomels Exp $	 */
