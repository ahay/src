#include <math.h>

#include <rsf.h>

#include "freqfilt.h"

static int nfft, nw;
static float complex *cdata;
static float *shape, *tmp;
static kiss_fftr_cfg forw, invs;

void freqfilt_init(int nfft1, int nw1)
{
    nfft = nfft1;
    nw = nw1;

    cdata = sf_complexalloc(nw);
    tmp = sf_floatalloc(nfft);
    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation problem",__FILE__);
}

void freqfilt_set(float *filt)
{
    shape = filt;
}

void freqfilt_close(void) {
    free(cdata);
    free(tmp);
    free(forw);
    free(invs);
}

void freqfilt_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
{
    int iw;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (iw=0; iw < nx; iw++) {
	tmp[iw] = adj? y[iw] : x[iw];
    }
    for (iw=nx; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    kiss_fftr(forw, tmp, (kiss_fft_cpx *) cdata);
    for (iw=0; iw < nw; iw++) {
	cdata[iw] *= shape[iw];
    }
    kiss_fftri(invs,(const kiss_fft_cpx *) cdata, tmp);

    for (iw=0; iw < nx; iw++) {	    
	if (adj) {
	    x[iw] += tmp[iw];
	} else {
	    y[iw] += tmp[iw];
	}
    } 
}

/* 	$Id$	 */
