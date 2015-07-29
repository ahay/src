/* Kolmogoroff spectral factorization */
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

#include <rsf.h>

#include "kolmog.h"

static kiss_fftr_cfg forw, invs;
static int nfft, nw, lag, shift;
static sf_complex *fft;

void kolmog_init(int n1, int lag1, int shift1)
/*< initialize with data length >*/
{
    nfft = n1;
    nw = nfft/2+1;
    lag = lag1;
    shift = shift1;

    fft = sf_complexalloc(nw);   

    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");
}

void kolmog_close(void)
/*< free allocated storage >*/
{
    free(fft);
    free(forw);
    free(invs);
}

void kolmog(float *trace)
/*< convert trace to minimum-phase >*/
{
    int i1;

    /* Fourier transform */
    kiss_fftr(forw,trace, (kiss_fft_cpx *) fft);
    for (i1=0; i1 < nw; i1++) {
	trace[i1] = 
	    crealf(fft[i1])*crealf(fft[i1])+
	    cimagf(fft[i1])*cimagf(fft[i1]);
    }

    kolmog2(trace);
}

void kolmog2(float *trace)
/*< convert Fourier-domain auto-correlation to minimum-phase >*/ 
{
    int i1;
    float asym;
    const double eps=1.e-32;

    for (i1=0; i1 < nw; i1++) {
	fft[i1] = sf_cmplx(log(trace[i1]+eps)/nfft,0.);
    }

    /* Inverse transform */
    kiss_fftri(invs,(const kiss_fft_cpx *) fft, trace);

    trace[0] *= 0.5;
    trace[nfft/2] *= 0.5;
    for (i1=1+nfft/2; i1 < nfft; i1++) {
	trace[i1] = 0.;
    }

    for (i1=1; i1 < lag; i1++) {
	asym = cosf(0.5*SF_PI*i1/(lag-1.0)); /* tapering weight */
	asym *= 0.5*trace[i1]*asym;
	trace[i1]      -= asym;
	trace[nfft-i1] += asym;
    }

    /* Fourier transform */
    kiss_fftr(forw,trace, (kiss_fft_cpx *) fft);
    
    for (i1=0; i1 < nw; i1++) {
#ifdef SF_HAS_COMPLEX_H
	fft[i1] = cexpf(fft[i1])/nfft;
#else
	fft[i1] = sf_crmul(cexpf(fft[i1]),1./nfft);
#endif
    }

    if (0 != shift) {
	for (i1=0; i1 < nw; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    fft[i1] *= cexpf(sf_cmplx(0.0,2*SF_PI*i1*shift/nfft));
#else
	    fft[i1] = sf_cmul(fft[i1],cexpf(sf_cmplx(0.0,2*SF_PI*i1*shift/nfft)));
#endif
	}
    }

    /* Inverse transform */
    kiss_fftri(invs,(const kiss_fft_cpx *) fft, trace);
}

