#include <rsf.h>

#include "kolmog.h"

static kiss_fftr_cfg forw, invs;
static int nfft, nw;
static float complex *fft;

void kolmog_init(int n1)
{
    nfft = n1;
    nw = nfft/2+1;
    fft = sf_complexalloc(nw);   

    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");
}

void kolmog_close(void)
{
    free(fft);
    free(forw);
    free(invs);
}

void kolmog(float *trace)
{
    int i1;

    /* Fourier transform */
    kiss_fftr(forw,trace, (kiss_fft_cpx *) fft);
    for (i1=0; i1 < nw; i1++) {
	trace[i1] = crealf(fft[i1]*conjf(fft[i1]));
    }

    kolmog2(trace);
}

void kolmog2(float *trace)
{
    int i1;
    const double eps=1.e-32;

    for (i1=0; i1 < nw; i1++) {
	fft[i1] = clog(trace[i1]+eps)/nfft;
    }

    /* Inverse transform */
    kiss_fftri(invs,(const kiss_fft_cpx *) fft, trace);

    trace[0] *= 0.5;
    trace[nfft/2] *= 0.5;
    for (i1=1+nfft/2; i1 < nfft; i1++) {
	trace[i1] = 0.;
    }

    /* Fourier transform */
    kiss_fftr(forw,trace, (kiss_fft_cpx *) fft);
    
    for (i1=0; i1 < nw; i1++) {
	fft[i1] = cexp(fft[i1])/nfft;
    }

    /* Inverse transform */
    kiss_fftri(invs,(const kiss_fft_cpx *) fft, trace);
}

