#include <float.h>
#include <math.h>

#include <rsf.h>

#include "xcorr.h"

float xcorr (int n1, int n2, const float *x1, const float *x2, 
	     int nc, float *xc)
{
    int is, ix, iy, imax, i2;
    float xmax, xp, xm, num, den, ds, shift, a1, a2;

    for (is = 0; is < nc; is++) {
	xc[is] = 0.; 
	a1 = 0.;
	a2 = 0.;
	for (i2=0; i2 < n2; i2++) {
	    for (ix=0; ix < n1; ix++) {
		iy = ix + (nc+1)/2 - is -1;
		if (iy >= 0 && iy < n1) {
		    xc[is] += x1[ix+i2*n1]*x2[iy+i2*n1];
		    a1 += x1[ix+i2*n1]*x1[ix+i2*n1];
		    a2 += x2[iy+i2*n1]*x2[iy+i2*n1];
		}
	    }
	}
	/* normalize */
	xc[is] /= sqrtf(a1*a2+FLT_EPSILON);
    }

    imax=(nc+1)/2;
    xmax=xc[imax];
    for (is = 0; is < nc; is++) {
	if (xc[is] > xmax) {
	    imax = is;
	    xmax = xc[is];
	}
    }
    
    shift = imax - (nc+1)/2;

    /* quadratic interpolation for sub-sample accuracy */
    if (imax > 0 && imax < nc-1) {
	xp = xc[imax+1];
	xm = xc[imax-1];
	num = 0.5*(xp-xm);
	den = xm+xp-2.*xmax;
	ds =  num*den/(den*den+FLT_EPSILON);
	if (fabsf(ds) < 1.) shift -= ds;
    }

    return shift;
}

/* 	$Id$	 */

