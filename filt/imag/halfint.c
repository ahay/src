#include "halfint.h"

#include <math.h>
#include <rsf.h>

static int nn, nw;
static float complex *cx, *cf;
static kiss_fftr_cfg forw, invs;

void halfint_init (bool adj, bool inv, int n, float rho)
{
    int i;
    float om;
    float complex cz;

    nn = n;
    nw = nn/2+1;

    cx = sf_complexalloc(nw);
    cf = sf_complexalloc(nw);

    forw = kiss_fftr_alloc(nn,0,NULL,NULL);
    invs = kiss_fftr_alloc(nn,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");

    for (i=0; i < nw; i++) {
	om = 2.*SF_PI*i/nn;
        if (!adj) om = - om;

	cz = cexpf(I*om);
	if (inv) {
	    cf[i] = csqrtf(1.-rho*cz);
	} else {
	    cf[i] = csqrtf(0.5*(1.+rho*cz)/(1.-rho*cz));
	}
    }
}

void halfint (float* x)
{
    int i;

    kiss_fftr(forw,x, (kiss_fft_cpx *) cx);
    for (i=0; i < nw; i++) {
	cx[i] *= cf[i]/nn;
    }
    kiss_fftri(invs,(const kiss_fft_cpx *) cx, x);
}

void halfint_close(void)
{
    free (cx);
    free (cf);
    free (forw);
    free (invs);
}

