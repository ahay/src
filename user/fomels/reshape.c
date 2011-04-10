/* Changing data spectrum by frequency-domain Gaussian scaling and smoothing */
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

#include "reshape.h"

static int nfft, nw;
static float dw;
static kiss_fft_cpx *cdata;
static kiss_fftr_cfg forw, invs;

void reshape_init(int n1, float d1 /* data length and sampling */)
/*< initialize >*/
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = n1;
    nw = nfft/2+1;
    dw = 1./(nfft*d1);

    cdata = (kiss_fft_cpx*) sf_complexalloc(nw);
    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);   
}

void reshape_close(void)
/*< free allocated storage >*/
{
    free (cdata);
    free (forw);
    free (invs);
}

void reshape (float m1, float a1 /* old parameters */, 
	      float m2, float a2 /* new parameters */, 
	      float* data        /* data - modified in place */) 
/*< smooth >*/
{
    int iw;
    float f;

    if (m1 < m2) return;

    kiss_fftr(forw,data,cdata);
    for (iw=0; iw < nw; iw++) {
	f = iw*dw;
	f *= f;
	f = (a2*m1)/(a1*m2)*expf(f*(1./m1-1./m2))/nfft;
	cdata[iw].r *= f;
	cdata[iw].i *= f;
    }
    kiss_fftri(invs,cdata,data);
}

/* 	$Id$	 */
