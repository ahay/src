#include <rsf.h>

#include "butter.h"

static float **den, *yy, sinw, cosw;
static int nn;

void butter_init(int npoly, int nx)
{
    nn = npoly;
    den = sf_floatalloc2(3,nn/2);
}

void butter_close(void)
{
    free (den);
}

void butter_set(float cutoff)
{
    int i;
    float arg, f;

    arg = 2.*SF_PI*cutoff;
    sinw = sinf(arg);
    cosw = cosf(arg);

    /*
    if (nn%2) {
	tanw = low? (1.+cosw)/sinw : sinw/(1.+cosw);
	den[0] = 1.+tanw;
	den[1] = 1.-tanw;
	n1 = 2;
    } else {
	den[0] = 1.;
	n1 = 1;
    }
    */

    for (i=0; i < nn/2; i++) {
	f = sinf(SF_PI*(2*i+1)/(2*nn))*sinw;
	den[i][0] = (1.+f);
	den[i][1] = -2.*cosw/(1.+f);
	den[i][2] = (1.-f)/(1.+f);
    }

    sinw = sinf(0.5*arg);
    cosw = cosf(0.5*arg);

    /*
      scale = low? sinf(0.5*arg): cosf(0.5*arg);
      scale = 1./(scale*scale);
    */
}

void butter (bool low, int ny, float *yy)
{
    /* convolution */
    for (iy=0; iy < ny; iy++) {
	yy[iy] = xx[iy];
    }
    for (ia=1; ia < na; ia++) { 
        for (iy=0; iy < ny; iy++) {
	    yy[iy+ia] += xx[iy] * num[ia];
	}
    }
    /* polynomial division */
    for (iy=0; iy < na-1; iy++) { /* lead-in terms */
	for (ia=1; ia <= iy; ia++) { 
	    yy[iy] -= den[ia] * yy[iy-ia];
	}
	yy[iy] /= den[0];
    }
    for (iy=na-1; iy < ny; iy++) { /* steady state */
	for (ia=1; ia < na; ia++) {
	    yy[iy] -= den[ia] * yy[iy-ia];
	}
	yy[iy] /= den[0];
    }
}
