#ifndef __SEAM_VERIFY
#define __SEAM_VERIFY

#include "utils.h"
#include "lci.h"
#include "su.h"
#include "header.h"
#include "segy.h"

/* trace comparison. computes nominal amplitude curve based on
   approximate reflection amplitude as function of source and receiver
   positions and time, returns relative max norm, normalized by this
   nominal amplitude series. The nominal reflection coefficient used
   to compute this reference amplitude is the fourth argument, rc.

   The direct wave amplitude, is amp*rdist/dist, where rdist is a
   reference distance - thus amp is meant to be the amplitude of the
   direct wave at the reference distance. The direct wave amplitude,
   at a distance equal to the two-way time multiplied by reference
   velocity, scaled by the refl coeff, is the proxy for reflected wave
   reference amplitude.

   To accomodate possibly different time steps and numbers of samples,
   use local cubic interpolation.
*/

int verifytrace(segy * trial, /* input trial trace */ 
		segy * comp,  /* comparison trace */
		int dim,      /* dimension of wave propagation */
		float tol,    /* tolerance for success */
		float rc,     /* nominal reflection coefficient */
		float cref,   /* nominal velocity */
		float amp,    /* nominal amplitude at reference distance */
                float rdist,  /* nominal reference distance for amplitue normalization */
		float * e);   /* scaled max error */

#endif
