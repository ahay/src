#include <math.h>

#include <rsf.h>

#include "window1.h"

static float h,dw;

void window1_init (float h_in, float dw_in)
{
    h = h_in;
    dw = dw_in;
}

int window1_apply (int iw, int w, const float* dat, 
		   bool left, bool right, float *win)
{
    int i, i0;
    float gain;

    i0 = 0.5+iw*dw;
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

/* 	$Id$	 */
