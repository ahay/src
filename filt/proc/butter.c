#include <rsf.h>

#include "butter.h"

void butter_set(bool low, float cutoff, int na, float *num, float *den)
{
    int i, j, n1, nn;
    float arg, ss, sinw, cosw, fact, scale;

    arg = 2.*SF_PI*cutoff;
    sinw = sinf(arg);
    cosw = cosf(arg);

    nn = na-1;

    if (nn%2) {
	if (low) {
	    fact = (1.+cosw)/sinw;
	    den[0] = 1.+fact;
	    den[1] = 1.-fact;
	} else {
	    fact = sinw/(1.+cosw);
	    den[0] = fact+1.;
	    den[1] = fact-1.;
	}
	n1 = 2;
    } else {
	den[0] = 1.;
	n1 = 1;
    }

    fact = low? sinf(0.5*arg): cosf(0.5*arg);
    fact = 1./(fact*fact);
    scale = 1.;

    for (i=0; i < nn/2; i++) {
	ss = sinf(SF_PI*(2*i+1)/(2*nn))*sinw;

	for (j=0; j < n1+2; j++) {
	    num[j] = 0.;
	}
	for (j=0; j < n1; j++) {
	    num[j]   += (1.+ss)*den[j];
	    num[j+1] -= 2.*cosw*den[j];
	    num[j+2] += (1.-ss)*den[j];
	}
	n1 += 2;
	for (j=0; j < n1; j++) {
	    den[j] = num[j];
	}
	scale *= fact;
    }
    
    for (j=0; j < na; j++) {
	den[j] *= scale;
    }
    
    num[0] = 1.;
    for (i=0; i < nn; i++) {
	num[i+1] = 1.;
	for (j=i-1; j >= 0; j--) {
	    num[j+1] += num[j];
	}
    }
    if (!low) {
	for (j=1; j < na; j+=2) {
	    num[j] = -num[j];
	}
    }
}

void butter (bool adj, int na, float *num, float *den, 
	     int nx, int ny, float *xx, float *yy)
{
    int ix, iy, ia;

    if (adj) {
	/* polynomial division */
	for (ix=ny-1; ix >=0; ix--) {
	    for (ia=1; ia < na; ia++) {
		iy = ix+ia;
		if (iy >= ny) break; 
		yy[ix] -= den[ia] * yy[iy];
	    }
	    yy[ix] /= den[0];
	}
	/* convolution */
	for (iy=0; iy < ny; iy++) {
	    xx[iy] = yy[iy];
	}
	for (ia=1; ia < na; ia++) { 
	    for (iy=ia; iy < ny; iy++) {
		ix=iy-ia;
		xx[ix] += yy[iy] * num[ia];
	    }
	}
    } else {
	/* convolution */
	for (ix=0; ix < nx; ix++) {
	    yy[ix] = xx[ix];
	}
	for (ia=1; ia < na; ia++) { 
	    for (iy=ia; iy < ny; iy++) {
		ix=iy-ia;
		yy[iy] += xx[ix] * num[ia];
	    }
	}
	/* polynomial division */
	for (iy=0; iy < ny; iy++) { /* lead-in terms */
	    for (ia=1; ia <= na; ia++) {
		ix = iy-ia;
		if (ix < 0) break;
		yy[iy] -= den[ia] * yy[ix];
	    }
	    yy[iy] /= den[0];
	}
    }
}
