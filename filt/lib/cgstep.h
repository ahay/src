#ifndef _sf_cgstep_h
#define _sf_cgstep_h

#include "c99.h"

/* cgstep
   --------
   A step of Claerbout's  conjugate-gradient iteration.
   nx - model size
   ny - data size
   x[nx] - current model
   g[nx] - gradient
   rr[ny] - residula (d - A x)
   gg[ny] - conjugate gradient */
void sf_cgstep( bool forget, int nx, int ny, 
		float* x, const float* g, float* rr, const float* gg);
void sf_cgstep_close (void);

#endif

/* 	$Id: cgstep.h,v 1.2 2004/03/13 06:00:15 fomels Exp $	 */
