#include <math.h>

#include <rsf.h>

#include "reshape.h"

static int nfft, nw;
static float dw;
static float complex *cdata;

int reshape_init(int n1, float d1)
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = sf_npfar(n1);
    nw = nfft/2+1;
    dw = 1./(nfft*d1);

    cdata = sf_complexalloc(nw);

    return nfft;
}

void reshape (float m1, float a1, float m2, float a2, float* data) 
{
    int iw;
    float f;

    if (m1 < m2) return;

    sf_pfarc (1,nfft,data,cdata);
    for (iw=0; iw < nw; iw++) {
	f = iw*dw;
	f *= f;
	cdata[iw] *= (a2*m1)/(a1*m2)*expf(f*(1./m1-1./m2))/nfft;
    }
    sf_pfacr (-1,nfft,cdata,data);
}
