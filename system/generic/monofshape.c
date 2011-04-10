/* Frequency-domain 1-D Gaussian shaper */
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
/*^*/

#include "monofshape.h"
#include "monof.h"

static int nfft, nw;
static kiss_fft_cpx *cdata;
static float *shape, *tmp, dw;
static kiss_fftr_cfg forw, invs;

void monofshape_init(int n1)
/*< initialize with data length >*/
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    cdata = (kiss_fft_cpx*) sf_complexalloc(nw);
    shape = sf_floatalloc(nw);
    tmp = sf_floatalloc(nfft);
    
    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation problem",__FILE__);
}

void monofshape_set(float a0       /* initial value for Gaussian */, 
		    int n          /* pattern length */, 
		    float* pattern /* data pattern [n] */, 
		    int niter      /* number of iterations */)
/*< set the shape >*/
{
    int iw, i0;
    float f, w, max, scale;
    
    for (iw=0; iw < n; iw++) {
	tmp[iw] = pattern[iw];
    }
    for (iw=n; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    scale = sqrtf(1./nfft); /* FFT scaling */ 

    kiss_fftr(forw, tmp,cdata);
    max = 0.;
    i0 = 0;
    for (iw=0; iw < nw; iw++) {
	f = sf_cabsf(cdata[iw])*scale;
	if (f > max) {
	    max = f;
	    i0 = iw;
	}
	tmp[iw] = f;
    }

    sf_warning("i0=%d max=%g",i0,max);

    f = monof(tmp,i0,niter,a0,nw,dw,true);

    for (iw=0; iw < nw; iw++) {
	w = (iw-i0)*dw;
	w *= w;
	shape[iw] = expf(-0.5*f*w)/nfft;
    } 
    
}

void monofshape_close(void) 
/*< free allocated storage >*/
{
    free(cdata);
    free(shape);
    free(tmp);
    free(forw);
    free(invs);
}

void monofshape_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator for shaping >*/
{
    int iw;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (iw=0; iw < nx; iw++) {	    
	tmp[iw] = adj? y[iw] : x[iw];
    }

    for (iw = nx; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    kiss_fftr(forw, tmp, cdata);
    for (iw=0; iw < nw; iw++) {
	cdata[iw] = sf_crmul(cdata[iw],shape[iw]);
    }
    kiss_fftri(invs, cdata, tmp);

    for (iw = 0; iw < nx; iw++) {
	if (adj) {
	    x[iw] += tmp[iw];
	} else {
	    y[iw] += tmp[iw];
	}
    }
}

/* 	$Id$	 */
