#include <rsf.h>

#include "kolmog.h"

void kolmog(int nfft, float *trace)
{
    int nw;
    float complex *fft;

    nw = nfft/2+1;
    fft = sf_complexalloc(nw);

    /* Fourier transform */
    sf_pfarc (1,nfft,trace,fft);

    kolmog2(nfft, nw, trace, fft);
    free(fft);
}

void kolmog2(int nfft, int nw, float *trace, float complex *fft)
{
    int i1;
    const double eps=1.e-32;

    for (i1=0; i1 < nw; i1++) {
	fft[i1] = clog(fft[i1]*conj(fft[i1])+eps)/nfft;
	sf_warning("got %d: (%g,%g)",i1,crealf(fft[i1]),cimagf(fft[i1]));
    }

    /* Inverse transform */
    sf_pfacr(-1,nfft,fft,trace);

    trace[0] *= 0.5;
    trace[nfft/2] *= 0.5;
    for (i1=1+nfft/2; i1 < nfft; i1++) {
	trace[i1] = 0.;
    }

    /* Fourier transform */
    sf_pfarc (1,nfft,trace,fft);

    for (i1=0; i1 < nw; i1++) {
	fft[i1] = cexp(fft[i1])/nfft;
	sf_warning("got %d: (%g,%g)",i1,crealf(fft[i1]),cimagf(fft[i1]));
    }

    /* Inverse transform */
    sf_pfacr(-1,nfft,fft,trace);
}



