#include <math.h>

#include <rsf.h>

#include "ricker.h"

static float complex *shape;

void ricker_init(int nfft   /* time samples */, 
		 float freq /* frequency */,
		 int order  /* derivative order */)
/*< initialize >*/
{
    int iw, nw;
    float dw, w;
    float complex cw;

    /* determine frequency sampling (for real to complex FFT) */
    nw = nfft/2+1;
    dw = 1./(nfft*freq);
 
    shape = sf_complexalloc(nw);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w *= w;

	switch (order) {
	    case 2: /* half-order derivative */
		cw = csqrtf((1.+I*iw)*2*SF_PI/nfft);
		shape[iw] = cw*w*expf(1-w)/nfft;
		break;
	    case 0:
	    default:
		shape[iw] = w*expf(1-w)/nfft;
		break;
	}
    }

    sf_freqfilt_init(nfft,nw);
    sf_freqfilt_cset(shape);
}

void ricker_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    sf_freqfilt_close();
}

/* 	$Id: ricker.c 694 2004-07-06 21:04:46Z fomels $	 */
