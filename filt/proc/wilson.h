#ifndef _wilson_h
#define _wilson_h

#include <rsf.h>

#include "helix.h"

/* Wilson-Burg spectral factorization */
void wilson_init( int nmax);
float wilson_factor(int niter, float s0, filter ss, 
		    filter aa, bool verb, float tol);
void wilson_close( void);

#endif
