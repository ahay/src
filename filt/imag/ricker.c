#include <math.h>

#include <rsf.h>

#include "ricker.h"

static float *shape;

void ricker_init(int nfft /* time samples */, 
		 float freq /* frequency */)
/*< initialize >*/
{
    int iw, nw;
    float dw, w, a;

    /* determine frequency sampling (for real to complex FFT) */
    nw = nfft/2+1;
    dw = 1./(nfft*freq);
    a = dw/powf(SF_PI,0.375);

    shape = sf_floatalloc(nw);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w *= w;
	shape[iw] = a*w*expf(-w);
    }

    sf_freqfilt_init(nfft,nw);
    sf_freqfilt_set(shape);
}

void ricker_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    sf_freqfilt_close();
}

/* 	$Id: ricker.c 694 2004-07-06 21:04:46Z fomels $	 */
