#include <rsf.h>

#include "butter.h"

struct Butter {
    bool low;
    int nn;
    float **den, mid;
};


butter butter_init(bool low, float cutoff, int nn)
{
    int j;
    float arg, ss, sinw, cosw, fact;
    butter bw;

    arg = 2.*SF_PI*cutoff;
    sinw = sinf(arg);
    cosw = cosf(arg);

    bw = (butter) sf_alloc (1,sizeof(*bw));
    bw->nn = nn;
    bw->low = low;
    bw->den = sf_floatalloc2(2,(nn+1)/2);

    if (nn%2) {
	if (low) {
	    fact = (1.+cosw)/sinw;
	    bw->den[nn/2][0] = 1./(1.+fact);
	    bw->den[nn/2][1] = 1.-fact;
	} else {
	    fact = sinw/(1.+cosw);
	    bw->den[nn/2][0] = 1./(fact+1.);
	    bw->den[nn/2][1] = fact-1.;
	}
    }

    fact = low? sinf(0.5*arg): cosf(0.5*arg);
    fact *= fact;
    
    for (j=0; j < nn/2; j++) {
	ss = sinf(SF_PI*(2*j+1)/(2*nn))*sinw;
	bw->den[j][0] = fact/(1.+ss);
	bw->den[j][1] = (1.-ss)/fact;
    }
    bw->mid = -2.*cosw/fact;

    return bw;
}

void butter_close(butter bw)
{
    free(bw->den[0]);
    free(bw->den);
    free(bw);
}

void butter_apply (const butter bw, int nx, float *x)
{
    int ix, j, nn;
    float d0, d1, d2, x0, x1, x2, y0, y1, y2;

    d1 = bw->mid;
    nn = bw->nn;

    if (nn%2) {
	d0 = bw->den[nn/2][0];
	d2 = bw->den[nn/2][1];
	x0 = y1 = 0.;
	for (ix=0; ix < nx; ix++) { 
	    x1 = x0; x0 = x[ix];
	    y0 = (bw->low)? 
		(x0 + x1 - d2 * y1)*d0:
		(x0 - x1 - d2 * y1)*d0;
	    y1 = y0; x[ix] = y0;
	}
    }

    for (j=0; j < nn/2; j++) {
	d0 = bw->den[j][0];
	d2 = bw->den[j][1];
	x1 = x0 = y1 = y2 = 0.;
	for (ix=0; ix < nx; ix++) { 
	    x2 = x1; x1 = x0; x0 = x[ix];
	    y0 = (bw->low)? 
		(x0 + 2*x1 + x2 - d1 * y1 - d2 * y2)*d0:
		(x0 - 2*x1 + x2 - d1 * y1 - d2 * y2)*d0;
	    y2 = y1; y1 = y0; x[ix] = y0;
	}
    }
}

