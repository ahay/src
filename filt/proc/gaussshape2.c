#include <math.h>

#include <rsf.h>

#include "gaussshape2.h"
#include "monof2.h"
#include "freqfilt2.h"

static int nfft, nw, nk;
static float **shape, dw, dk, w0;

void gaussshape2_init(int n1, int n2)
{
    /* determine frequency sampling (for real to complex FFT) */
    nw = sf_npfa(n1);
    dw = 2.*SF_PI/nw;
    w0 = -SF_PI;

    /* determine wavenumber sampling (for complex FFT) */
    nfft = sf_npfar(n2);
    nk = nfft/2+1;
    dk = 2.*SF_PI/nfft;

    shape = sf_floatalloc2(nw,nk);
    freqfilt2_init(n1, n2, nfft, nw, nk);
    freqfilt2_set(shape);
}

void gaussshape2_set(float* a, const float* pattern, int niter)
{
    freqfilt2_spec(pattern,shape);
    monof2(shape, niter, a, nw, dw, w0, nk, dk, 0., true);
    gaussshape2_set2(a);
}

void gaussshape2_set2(const float* a) {
    int ik, iw;
    float w, w2, k, k2, wk;
    
    for (ik=0; ik < nk; ik++) {
	k = ik*dk;
	k2 = k*k;
	for (iw=0; iw < nw; iw++) {
	    w = w0+iw*dw;
	    w2 = w*w;
	    wk = w*k;
	    shape[ik][iw] = expf(-0.5*(a[0]*w2+a[1]*wk+a[2]*k2))/(nfft*nw);
	}
    }
}

void gaussshape2_close(void) {
    free(shape[0]);
    free(shape);
    freqfilt2_close();
}

/* 	$Id: gaussshape2.c,v 1.2 2004/02/26 05:16:08 fomels Exp $	 */
