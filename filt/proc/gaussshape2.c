#include <math.h>

#include <rsf.h>

#include "gaussshape2.h"
#include "monof2.h"
#include "freqfilt2.h"

static int nfft, nw, nk;
static float **shape, dw, dk, k0;

void gaussshape2_init(int n1, int n2)
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = n1;
    if (n1%2) nfft++;
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    /* determine wavenumber sampling (for complex FFT) */
    nk = n2;
    dk = 2.*SF_PI/n2;
    k0 = -SF_PI;

    shape = sf_floatalloc2(n2,nw);
    freqfilt2_init(n1, n2, nw);
    freqfilt2_set(shape);
}

void gaussshape2_set(float* a, const float* pattern, int niter)
{
    freqfilt2_spec(pattern,shape);
    monof2(shape, niter, a, nk, dk, k0, nw, dw, 0., true);
    gaussshape2_set2(a);
}

void gaussshape2_set2(const float* a) {
    int ik, iw;
    float w, w2, k, k2, wk;

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w2 = w*w;
	for (ik=0; ik < nk; ik++) {
	    k = k0+ik*dk;
	    k2 = k*k;
	    wk = w*k;
	    
	    shape[iw][ik] = expf(-0.5*(a[0]*w2+a[1]*wk+a[2]*k2))/(nfft*nk);
	}
    }
}

void gaussshape2_close(void) {
    free(shape[0]);
    free(shape);
    freqfilt2_close();
}

/* 	$Id$	 */
