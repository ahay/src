#include <math.h>

#include <rsf.h>

#include "window3.h"

static int *w, *nw, *n;
static float *h;

void window3_init (int* w_in, int* nw_in, int* n_in, float* h_in)
{
    w = w_in; 
    nw = nw_in; 
    n = n_in; 
    h = h_in;
}

void window3_apply (const int* i, float*** dat, 
		    bool near, bool far, 
		    bool left, bool right, 
		    bool top, bool bottom,
		    int* i0, float*** win)
{
    int i1, i2, i3, j;
    float gain1, gain2, gain3;

    for (j=0; j < 3; j++) {
	i0[j] = 0.5+(n[j]-w[j])*i[j]/(nw[j]-1.);
    }

    for (i3=0; i3 < w[2]; i3++) {
	if (near && i3+1 <= h[2]) {
	    gain3 = sinf(0.5*SF_PI*i3/h[2]);
	} else if (far && i3 >= w[2]-h[2]) {
	    gain3 = sinf(0.5*SF_PI*(w[2]-i3-1)/h[2]);
	} else {
	    gain3 = 1.;
	}
	for (i2=0; i2 < w[1]; i2++) {
	    if (left && i2+1 <= h[1]) {
		gain2 = sinf(0.5*SF_PI*i2/h[1]);
	    } else if (right && i2 >= w[1]-h[1]) {
		gain2 = sinf(0.5*SF_PI*(w[1]-i2-1)/h[1]);
	    } else {
		gain2 = 1.;
	    }
	    
	    gain2 *= gain3;

	    for (i1=0; i1 < w[0]; i1++) {
		if (top && i1+1 <= h[0]) {
		    gain1 = sinf(0.5*SF_PI*i1/h[0]);
		} else if (bottom && i1 >= w[0]-h[0]) {
		    gain1 = sinf(0.5*SF_PI*(w[0]-i1-1)/h[0]);
		} else {
		    gain1 = 1.;
		}
		
		gain1 *= gain2;
		gain1 *= gain1;
          
		win[i3][i2][i1] = gain1*dat[i0[2]+i3][i0[1]+i2][i0[0]+i1];
	    }
	}
    }
}

/* 	$Id: window3.c,v 1.2 2004/06/16 17:55:15 fomels Exp $	 */
