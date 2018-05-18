/* Frequency-domain filtering. */
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

#include "freqfilt.h"

#include "_bool.h"
#include "kiss_fft.h"
/*^*/

#include "alloc.h"
#include "komplex.h"
#include "error.h"
#include "adjnull.h"
#include "kiss_fftr.h"
#include "_kiss_fft_guts.h"

static int nfft, nw;
static kiss_fft_cpx *cdata, *shape=NULL;
static float *tmp;
static kiss_fftr_cfg forw, invs;

void sf_freqfilt_init(int nfft1 /* time samples (possibly padded) */, 
		      int nw1   /* frequency samples */)
/*< Initialize >*/
{
    nfft = nfft1;
    nw = nw1;

    cdata = (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
    tmp = sf_floatalloc(nfft);
    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation problem",__FILE__);
}

void sf_freqfilt_set(float *filt /* frequency filter [nw] */)
/*< Initialize filter (zero-phase) >*/
{
    int iw;
    
    if (NULL==shape) shape = (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));

    for (iw=0; iw < nw; iw++) {
	shape[iw].r = filt[iw];
	shape[iw].i = 0.;
    }
}

/* #ifndef __cplusplus */
/*^*/

void sf_freqfilt_cset(kiss_fft_cpx *filt /* frequency filter [nw] */)
/*< Initialize filter >*/
{
    shape = filt;
}

void sf_freqfilt_close(void) 
/*< Free allocated storage >*/
{
    free(cdata);
    free(tmp);
    free(forw);
    free(invs);
}

void sf_freqfilt(int nx, float* x)
/*< Filtering in place >*/
{
    kiss_fft_cpx c;
    int iw;

    for (iw=0; iw < nx; iw++) {
	tmp[iw] = x[iw];
    }
    for (iw=nx; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    kiss_fftr(forw, tmp, cdata);
    for (iw=0; iw < nw; iw++) {
	C_MUL(c,cdata[iw],shape[iw]);
	cdata[iw]=c;
    }
    kiss_fftri(invs, cdata, tmp);

    for (iw=0; iw < nx; iw++) {
	x[iw] = tmp[iw];
    } 
}

void sf_freqfilt_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< Filtering as linear operator >*/
{
    kiss_fft_cpx c;
    int iw;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (iw=0; iw < nx; iw++) {
	tmp[iw] = adj? y[iw] : x[iw];
    }
    for (iw=nx; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    kiss_fftr(forw, tmp, cdata);
    for (iw=0; iw < nw; iw++) {
        if (adj) {
	    C_MUL(c,cdata[iw],sf_conjf(shape[iw]));
        } else {
            C_MUL(c,cdata[iw],shape[iw]);
        }
	cdata[iw]=c;
    }
    kiss_fftri(invs, cdata, tmp);

    for (iw=0; iw < nx; iw++) {	    
	if (adj) {
	    x[iw] += tmp[iw];
	} else {
	    y[iw] += tmp[iw];
	}
    } 
}

/* #endif */

/* 	$Id$	 */
