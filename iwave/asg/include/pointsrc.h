#ifndef __SEAM_POINT_SOURCE__
#define __SEAM_POINT_SOURCE__

#include "utils.h"
#include "su.h"
#include "header.h"
#include "segy.h"
#include "cubic.h"
#include "sgn.h"
#include "gauss.h"
#include "defaults.h"
#include "traceio.h"
#include "model.h"

/** Structure for dilatational point source terminator. Stores number
    of wavelet samples at simulation sample rate, start index for
    source action, source location information (cell index, position
    in cell), sampling order, and buoyancy values near source, along
    with the wavelet time series (at simulation sample rate). Wavelet
    is either 
       (i)  read from file, or
       (ii) extracted from target pulse, read from file, or
       (iii)Gaussian with zero, min, or max phase, producing
            a propagating Ricker pulse in 3D.
*/
typedef struct {
  int istart;          /* start index for source */
  IPNT is;             /* source indices */
  RPNT rs;             /* source in-cell offsets */
  ireal *w;            /* wavelet */
  int n;               /* wavelet length */
  int order;           /* sampling order */
  ireal scramp;	       /* scale factor for discrete delta_function */
  FILE *fpsrc;         /* stream for read from file */
  int idbg;            /* flag for optional output of working wavelet */
  FILE * fpdbg;	       /* stream for optional output of working wavelet */
  ireal bm;
} POINTSRC;

/* no need for default constructor, since there is no
   abstraction at this level */
int pointsrc_init(POINTSRC * t, 
		  IMODEL * m, 
		  PARARRAY * par, 
		  tracegeom * tg, 
		  FILE * stream);
int pointsrc_destroy(POINTSRC * t);
int pointsrc_run(POINTSRC * t, IMODEL * m);
void pointsrc_fprint(POINTSRC const * t, FILE * fp);

#endif
