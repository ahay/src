/* Trapezoidal filter */
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

#include "trapez.h"

static float *shape;

void trapez_init(int nfft   /* time samples */, 
		 float *freq /* frequencies */)
/*< initialize >*/
{
    int iw, nw;
    float dw, w, t;

    /* determine frequency sampling (for real to complex FFT) */
    nw = nfft/2+1;
    dw = 1./nfft;
 
    shape = sf_floatalloc(nw);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw;

	if (w < freq[0]) {
	    shape[iw] = 0.;
	} else if (w < freq[1]) {
	    t = sinf(0.5*SF_PI*(freq[1]-w)/(freq[1]-freq[0]));
	    shape[iw] = 1.-t*t;
	} else if (w < freq[2]) {
	    shape[iw] = 1.;
	} else if (w < freq[3]) {
	    t = sinf(0.5*SF_PI*(freq[3]-w)/(freq[3]-freq[2]));
	    shape[iw] = t*t;
	} else {
	    shape[iw] = 0.;
	}
    }

    sf_freqfilt_init(nfft,nw);
    sf_freqfilt_set(shape);
}

void trapez_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    sf_freqfilt_close();
}

/* 	$Id: trapez.c 694 2004-07-06 21:04:46Z fomels $	 */
