/* Ricker filter */
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

#include "ricker.h"

static kiss_fft_cpx *shape;

static kiss_fft_cpx my_csqrtf (float re, float im)
{
    float d, r, s;
    kiss_fft_cpx v;

    if (im == 0) {
      if (re < 0) {
	  v.r = 0.;
	  v.i = copysignf (sqrtf (-re), im);
      } else {
	  v.r =  fabsf (sqrtf (re));
	  v.i =  copysignf (0, im);
      }
    } else if (re == 0) {
	r = sqrtf (0.5 * fabsf (im));
	v.r = r;
	v.i = copysignf (r, im);
    } else {
	d = hypotf (re, im);
	/* Use the identity   2  Re res  Im res = Im x
	   to avoid cancellation error in  d +/- Re x.  */
	if (re > 0) {
	    r = sqrtf (0.5 * d + 0.5 * re);
	    s = (0.5 * im) / r;
        } else {
	    s = sqrtf (0.5 * d - 0.5 * re);
	    r = fabsf ((0.5 * im) / s);
        }
	v.r = r;
	v.i = copysignf (s, im);
    }
    return v;
}

void ricker_init(int nfft   /* time samples */, 
		 float freq /* frequency */,
		 int order  /* derivative order */)
/*< initialize >*/
{
    int iw, nw;
    float dw, w;
    kiss_fft_cpx cw;

    /* determine frequency sampling (for real to complex FFT) */
    nw = nfft/2+1;
    dw = 1./(nfft*freq);
 
    shape = (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;
	w *= w;

	switch (order) {
	    case 2: /* half-order derivative */
		cw = my_csqrtf(2*SF_PI/nfft,iw*2*SF_PI/nfft);
		shape[iw].r = cw.r*w*expf(1-w)/nfft;
		shape[iw].i = cw.i*w*expf(1-w)/nfft;
		break;
	    case 0:
	    default:
		shape[iw].r = w*expf(1-w)/nfft;
		shape[iw].i = 0.;
		break;
	}
    }

    sf_freqfilt_init(nfft,nw);
    sf_freqfilt_cset(shape);
}

void ricker_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    sf_freqfilt_close();
}

/* 	$Id: ricker.c 694 2004-07-06 21:04:46Z fomels $	 */
