#include <math.h>

#include <rsf.h>

#include "window1.h"

static int w,nw,n;
static float h;

void window1_init (int w_in, int nw_in, int n_in, float h_in)
{
    w = w_in; 
    nw = nw_in; 
    n = n_in; 
    h = h_in;
}

int window1_apply (int iw, float* dat, bool left, bool right, float *win)
{
    int i, i0;
    float gain;

    i0 = 0.5+(n-w)*iw/(nw-1.);
    for (i=0; i < w; i++) {
	if (left && i < h) {
	    gain = sinf(0.5*SF_PI*i/h);
	    win[i] = gain*gain*dat[i0+i];
	} else if (right && i >= w-h) {
	    gain =  sinf(0.5*SF_PI*(w-i-1)/h);
	    win[i] = gain*gain*dat[i0+i];
	} else {
	    win[i] = dat[i0+i];
	}
    }

    return i0;
}

