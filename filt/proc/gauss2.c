#include <math.h>

#include <rsf.h>

#include "gauss2.h"
#include "freqfilt2.h"

static float **shape;

void gauss2_init(int n1, int n2, float f1, float f2)
{
    int ik, iw, nfft, nw, nk;
    float dw, w, w2, dk, k, w0, k2;

    /* determine frequency sampling (for real to complex FFT) */
    nw = sf_npfa(n1);
    dw = 2.*SF_PI/nw;
    w0 = -SF_PI;

    /* determine wavenumber sampling (for complex FFT) */
    nfft = sf_npfar(n2);
    nk = nfft/2+1;
    dk = 2.*SF_PI/nfft;

    shape = sf_floatalloc2(nw,nk);

    f1 = sqrtf((f1*f1-1.)/12.);
    f2 = sqrtf((f2*f2-1.)/12.);

    for (ik=0; ik < nk; ik++) {
	k = ik*dk;
	k2 = f2*k;
	k2 *= k2;
	for (iw=0; iw < nw; iw++) {
	    w = w0+iw*dw;
	    w2 = f1*w;
	    w2 *= w2;
	    shape[ik][iw] = expf(-k2-w2)/(nfft*nw);
	}
    }

    freqfilt2_init(n1,n2,nfft,nw,nk);
    freqfilt2_set(shape);
}

void gauss2_close(void) {
    free(shape[0]);
    free(shape);
    freqfilt2_close();
}

/* 	$Id: gauss2.c,v 1.2 2004/02/26 05:16:08 fomels Exp $	 */
