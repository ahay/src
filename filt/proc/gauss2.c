#include <math.h>

#include <rsf.h>

#include "gauss2.h"
#include "freqfilt2.h"

static float **shape;

void gauss2_init(int n1, int n2, float f1, float f2)
{
    int ik, iw, nfft, nw;
    float dw, w, w2, dk, k, k0, k2;

    /* determine frequency sampling (for real to complex FFT) */
    nfft = n1;
    if (n1%2) nfft++;
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;
    
    /* determine wavenumber sampling (for complex FFT) */
    dk = 2.*SF_PI/n2;
    k0 = -SF_PI;

    shape = sf_floatalloc2(n2,nw);

    f1 = sqrtf((f1*f1-1.)/12.);
    f2 = sqrtf((f2*f2-1.)/12.);
 
    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w2 = f1*w;
	w2 *= w2;
	for (ik=0; ik < n2; ik++) {
	    k = k0+ik*dk;
	    k2 = f2*k;
	    k2 *= k2;
	    shape[iw][ik] = expf(-k2-w2)/(nfft*n2);
	}
    }

    freqfilt2_init(n1,n2,nw);
    freqfilt2_set(shape);
}

void gauss2_close(void) {
    free(shape[0]);
    free(shape);
    freqfilt2_close();
}

/* 	$Id$	 */
