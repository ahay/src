#ifndef _sf_ccgstep_h
#define _sf_ccgstep_h

#include "c99.h"

/* ccgstep
   --------
   A step of Claerbout's  conjugate-gradient iteration for complex operators.
   nx - model size
   ny - data size
   x[nx] - current model
   g[nx] - gradient
   rr[ny] - residula (d - A x)
   gg[ny] - conjugate gradient */
void sf_ccgstep( bool forget, int nx, int ny, 
		float complex* x, const float complex* g, 
		float complex* rr, const float complex* gg);
void sf_ccgstep_close (void);

#endif

/* 	$Id: ccgstep.h,v 1.1 2004/03/13 06:11:03 fomels Exp $	 */
