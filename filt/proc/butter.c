#include <rsf.h>

#include "butter.h"
#include "kolmog.h"

static const int nw, na;
static float *cw;

void butter_init(int nw_in)
{
    nw=sf_npfar(nw_in);
    cw = sf_floatalloc(nw);
    kolmog_init(nw);
}

void butter_close(void)
{
    free(cw);
    kolmog_close();
}

void butter_set(bool low, float cutoff, int npoly, float *num, float *den)
{
    int nn, i;
    real arg, tancut;

    na = npoly;
    nn = npoly - 1;
    tancut = 1./tanf(cutoff*SF_PI);
    for (i=0; i < nw; i++) {
        arg = SF_PI * i / nw;
	cx[i] = powf(cosf(arg),2*nn) + powf(sinf(arg) * tancut,2*nn);
    }
    kolmog(nw, cx); /* spectral factorization */
    for (i=0; i < npoly; i++) {
	den[i] = cx[i]; /* denominator */
    }
    for (i=0; i < nw; i++) {
        arg = SF_PI * i / nw;
	cx[i] = low? powf(cosf(arg),2*nn): powf(sinf(arg) * tancut,2*nn);
    }
    kolmog(nw, cx); /* spectral factorization */
    for (i=0; i < npoly; i++) {
	num[i] = cx[i]; /* numerator */
    }
}

void butter (int nx, int ny, const float *num, const float *den, 
	     const float *xx, float *yy)
{
    int ia, iy;

    /* convolution */
    for (iy=0; iy < ny; iy++) {
	yy[iy] = 0.;
    }
    for (ia=0; ia < na; ia++) { 
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
