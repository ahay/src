/* Butterworth filtering. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* Implementation is inspired by D. Hale and J.F. Claerbout, 1983, Butterworth
 * dip filters: Geophysics, 48, 1033-1038. */
    
#include "_bool.h"
/*^*/

#include "butter.h"
#include "_defs.h"
#include "c99.h"
#include "alloc.h"

#ifndef _sf_butter_h

typedef struct Sf_Butter *sf_butter;
/* abstract data type */
/*^*/

#endif

struct Sf_Butter {
    bool low;
    int nn;
    float **den, mid;
};


sf_butter sf_butter_init(bool low     /* low-pass (or high-pass) */, 
		   float cutoff /* cut off frequency */, 
		   int nn       /* number of poles */)
/*< initialize >*/
{
    int j;
    float arg, ss, sinw, cosw, fact;
    sf_butter bw;

    arg = 2.*SF_PI*cutoff;
    sinw = sinf(arg);
    cosw = cosf(arg);

    bw = (sf_butter) sf_alloc (1,sizeof(*bw));
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

void sf_butter_close(sf_butter bw)
/*< Free allocated storage >*/
{
    free(bw->den[0]);
    free(bw->den);
    free(bw);
}

void sf_butter_apply (const sf_butter bw, int nx, float *x /* data [nx] */)
/*< filter the data (in place) >*/
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
	    x[ix] = y1 = y0;
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
	    y2 = y1; x[ix] = y1 = y0;
	}
    }
}

void sf_reverse (int n1, float* trace)
/*< reverse a trace >*/
{
    int i1;
    float t;

    for (i1=0; i1 < n1/2; i1++) { 
        t=trace[i1];
        trace[i1]=trace[n1-1-i1];
        trace[n1-1-i1]=t;
    }
}

