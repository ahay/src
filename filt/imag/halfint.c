#include "halfint.h"

#include <math.h>
#include <rsf.h>

static int n, nn, nw;
static float *pad;
static float complex *cx, *cf;

void halfint_init (bool adj, bool inv, int n_in, float rho)
{
    int i;
    float om;
    float complex cz;

    n = n_in;

    nn = sf_npfar(n);
    nw = nn/2+1;

    pad = sf_floatalloc(nn);
    cx = sf_complexalloc(nw);
    cf = sf_complexalloc(nw);

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

    for (i=0; i < n; i++) {
	pad[i]=x[i];
    }
    for (i=n; i < nn; i++) {
	pad[i]=0.;
    }

    sf_pfarc(1,nn,pad,cx);
    for (i=0; i < nw; i++) {
	cx[i] *= cf[i]/nn;
    }
    sf_pfacr(-1,nn,cx,pad);

    for (i=0; i < n; i++) {
	x[i] = pad[i];
    }
}

void halfint_close(void)
{
    free (pad);
    free (cx);
    free (cf);
}

