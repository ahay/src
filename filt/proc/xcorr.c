#include <float.h>
#include <math.h>

#include <rsf.h>

#include "xcorr.h"

static int nc, nx, n2;
static float *xc;

void xcorr_init (int nx_in, int n2_in, int maxshift)
{
    nx = nx_in;
    n2 = n2_in;
    nc = 2*maxshift-1;
    xc = sf_floatalloc(nc);
}
  
void xcorr_close (void)
{
    free (xc);
}

float xcorr (const float *x1, const float *x2)
{
    int is, ix, iy, imax, i2;
    float xmax, xp, xm, num, den, ds, shift;

    for (is = 0; is < nc; is++) {
	xc[is] = 0.;
	for (i2=0; i2 < n2; i2++) {
	    for (ix=0; ix < nx; ix++) {
		iy = ix + (nc+1)/2 - is -1;
		if (iy >= 0 && iy < nx) 
		    xc[is] += x1[ix+i2*nx]*x2[iy+i2*nx];
	    }
	}
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
