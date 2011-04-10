/* Gaussian shaping in 2-D */
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

#include <math.h>

#include <rsf.h>

#include "gaussshape2.h"
#include "monof2.h"

static int nfft, nw, nk;
static float **shape, **shape2, dw, dk, k0;

void gaussshape2_init(int n1, int n2)
/*< initialize with data size >*/
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    /* determine wavenumber sampling (for complex FFT) */
    nk = n2;
    dk = 2.*SF_PI/n2;
    k0 = -SF_PI;

    shape  = sf_floatalloc2(n2,nw);
    shape2 = sf_floatalloc2(n2,nw);

    sf_freqfilt2_init(n1, n2, nw);
    sf_freqfilt2_set(shape);
}

void gaussshape2_set(float* a             /* shaper */, 
		     const float* pattern /* training data */, 
		     int niter            /* number of iterations */,
		     int nliter           /* number of reweighting iterations */)
/*< estimate shaping >*/
{
    sf_freqfilt2_spec(pattern,shape);
    monof2(shape,shape2,nliter,niter, a, nk, dk, k0, nw, dw, 0., true);
    gaussshape2_set2(a);
}

void gaussshape2_set2(const float* a) 
/*< set shaper >*/
{
    int ik, iw;
    float w, w2, k, k2, wk;

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w2 = w*w;
	for (ik=0; ik < nk; ik++) {
	    k = k0+ik*dk;
	    k2 = k*k;
	    wk = w*k;
	    
	    shape[iw][ik] = expf(-0.5*(a[0]*w2+a[1]*wk+a[2]*k2))/(nfft*nk);
	}
    }
}

void gaussshape2_close(void) 
/*< free allocated storage >*/
{
    free(shape[0]);
    free(shape);
    sf_freqfilt2_close();
}

/* 	$Id$	 */
