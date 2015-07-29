/* 1-D Gaussian frequency-domain smoothing */
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

#include "gauss.h"

static float *shape;

void gauss_init(int nx     /* trace length */, 
		float rect /* smoothing length */)
/*< initialize (call sf_freqfilt afterwards) >*/
{
    int iw, nw, nfft;
    float dw, w;

    /* determine frequency sampling (for real to complex FFT) */
    nfft = nx%2? nx+1:nx;
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    rect = sqrtf((rect*rect-1.)/12.);

    shape = sf_floatalloc(nw);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw*rect;
	w *= w;
	shape[iw] = expf(-w)/nfft;
    }

    sf_freqfilt_init(nfft,nw);
    sf_freqfilt_set(shape);
}

void gauss_close(void) 
/*< free allocated storage >*/
{
    free(shape);
    sf_freqfilt_close();
}

/* 	$Id: gauss.c 7107 2011-04-10 02:04:14Z ivlad $	 */
