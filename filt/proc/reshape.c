#include <math.h>

#include <rsf.h>

#include "reshape.h"

static int nfft, nw;
static float dw;
static float complex *cdata;
static kiss_fftr_cfg forw, invs;

void reshape_init(int n1, float d1)
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = n1;
    nw = nfft/2+1;
    dw = 1./(nfft*d1);

    cdata = sf_complexalloc(nw);
    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);   
}

void reshape_close(void)
{
    free (cdata);
    free (forw);
    free (invs);
}

void reshape (float m1, float a1, float m2, float a2, float* data) 
{
    int iw;
    float f;

    if (m1 < m2) return;

    kiss_fftr(forw,data, (kiss_fft_cpx *) cdata);
    for (iw=0; iw < nw; iw++) {
	f = iw*dw;
	f *= f;
	cdata[iw] *= (a2*m1)/(a1*m2)*expf(f*(1./m1-1./m2))/nfft;
    }
    kiss_fftri(invs,(const kiss_fft_cpx *) cdata, data);
}

/* 	$Id$	 */
