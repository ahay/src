#include <math.h>

#include <rsf.h>

#include "ricker.h"
#include "freqfilt.h"

static float *shape;

void ricker_init(int nfft /* time samples */, 
		 float freq /* frequency */)
/*< initialize >*/
{
    int iw, nw;
    float dw, w, a;

    /* determine frequency sampling (for real to complex FFT) */
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    shape = sf_floatalloc(nw);

    a = sqrtf(SF_PI)/(2*nfft*freq);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw*freq;
	w *= w;
	shape[iw] = a*w*expf(-w);
    }

    freqfilt_init(nfft,nw);
    freqfilt_set(shape);
}

void ricker_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    freqfilt_close();
}

/* 	$Id: ricker.c 694 2004-07-06 21:04:46Z fomels $	 */
