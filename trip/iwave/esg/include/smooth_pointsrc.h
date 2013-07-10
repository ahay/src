#ifndef __SEAM_SMOOTH_POINT_SOURCE__
#define __SEAM_SMOOTH_POINT_SOURCE__

#include "utils.h"
#include "su.h"
#include "header.h"
#include "segy.h"
#include "cubic.h"
#include "esgn.h"
#include "gauss.h"
#include "defaults.h"
#include "traceio.h"
#include "model.h"

/** 
    Structure for smooth point source. Stores number
    of wavelet samples at simulation sample rate, start index for
    source action, source location information (cell index, position
    in cell), sampling order, and bouyancy values near source, along
    with the wavelet time series (at simulation sample rate). Wavelet
    is either 
       (i)  read from file, or
       (ii) extracted from target pulse, read from file, or
       (iii)Gaussian with zero, min, or max phase, producing
            a propagating Ricker pulse in 3D.
            
    Source function for the pressure equation:
      -\nabla \phi(r) \cdot v0(r,t)
    Source function for the velocity equations:
      -\nabla \phi(r) p0(r,t),
    where  
      p0, v0, are the slns of the homogeneous acoustic problem given by
       eqns. (5) and (8) in the tech report by WWS&TV;
       
      phi = 1 for r < rad, phi = gauss(phipeak, rad) for r >= rad.
      
      Therefore,
      \nabla phi(r) = x phi'(r) / r and
      phi'(r) = - 2 * phipeak^2 * pi^2 * (r - rad) * exp(-phipeak^2 * pi^2 * (r - rad)^2 );
       
    
*/
typedef struct {
  int istart;          /* start index for source */
  RPNT xs;             /* physical coordinates of the source */
  ireal bou;           /* bouyancy at source location */
	ireal c;             /* sound velocity at source location */
	ireal fpeak;				 /* peak frequency for the ricker wavelet */ 
	IPNT ixs;						 /* beginning index for the cutoff function */
	IPNT ixe;						 /* end index for the cutoff function */
	FILE *fpsrc;         /* stream for read from file */
	ireal *w;            /* wavelet */
	ireal *w1;           /* integrated wavelet in option I(a), differentiated wavelet in option I(b) */
  int n;               /* wavelet length */
  ireal dt;            /* wavelet step   */
  ireal t0;            /* wavelet start  */
	int idbg;            /* flag for optional output of working wavelet */
  FILE * fpdbg;	       /* stream for optional output of working wavelet */
	int tmpflag;         /* 3 - option I(a), 2 - option I(b), 1 - option II (Ricker) */
	ireal scramp; 
  ireal phipeak;       /* half support for gaussian for the cutoff function */ 
  ireal rad;           /* radius of the flat area in the cutoff function    */
  int sn;              /* time (in iterations) that wavelet needs to get out 
                          of the homogeneous area around the source, computed as: 
                          rad + half support of the cutoff funcion + pulse width / sound velocity /dt
                        */
	
} SPOINTSRC;

/* no need for default constructor, since there is no
   abstraction at this level */
int spointsrc_init(SPOINTSRC * t, 
		  IMODEL * m, 
		  PARARRAY * par, 
		  tracegeom * tg, 
		  FILE * stream);
int spointsrc_destroy(SPOINTSRC * t);
int spointsrc_run(SPOINTSRC * t, IMODEL * m);
void spointsrc_fprint(SPOINTSRC const * t, FILE * fp);

#endif
