#ifndef _cgstep_h
#define _cgstep_h

#include <rsf.h>

/* cgstep
   --------
   A step of Claerbout's  conjugate-gradient interation.
   nx - model size
   ny - data size
   x[nx] - current model
   g[nx] - gradient
   rr[ny] - residula (d - A x)
   gg[ny] - conjugate gradient */
void cgstep( bool forget, int nx, int ny, 
	       float* x, const float* g, float* rr, const float* gg);
void cgstep_close (void);

#endif
