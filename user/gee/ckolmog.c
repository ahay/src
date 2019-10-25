/* complex Kolmogoroff spectral factorization */
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

#include "ckolmog.h"

static kiss_fft_cfg forw, invs;
static int nfft;
static sf_complex *fft;

void ckolmog_init(int n1)
/*< initialize with data length >*/
{
    nfft = n1;

    fft = sf_complexalloc(nfft);

    forw = kiss_fft_alloc(nfft,0,NULL,NULL);
    invs = kiss_fft_alloc(nfft,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");
}

void ckolmog_close(void)
/*< free allocated storage >*/
{
    free(fft);
    free(forw);
    free(invs);
}

void ckolmog(float *aut, sf_complex *trace)
/*< convert trace to minimum-phase >*/
{
    int i1;

    /* Fourier transform */
    kiss_fft(forw,(kiss_fft_cpx *) trace, (kiss_fft_cpx *) fft);
    for (i1=0; i1 < nfft; i1++) {
	aut[i1] = 
	    crealf(fft[i1])*crealf(fft[i1])+
	    cimagf(fft[i1])*cimagf(fft[i1]);
    }

    ckolmog2(aut, trace);
}

void ckolmog2(float *aut, sf_complex *trace)
/*< convert Fourier-domain auto-correlation to minimum-phase >*/ 
{
    int i1;
    const double eps=1.e-32;

    for (i1=0; i1 < nfft; i1++) {
	fft[i1] = sf_cmplx(log(aut[i1]+eps)/nfft,0.);
    }

    /* Inverse transform */
    kiss_fft(invs,(const kiss_fft_cpx *) fft, (kiss_fft_cpx *) trace);

    trace[0] *= 0.5;
    trace[nfft/2] *= 0.5;
    for (i1=1+nfft/2; i1 < nfft; i1++) {
	trace[i1] = 0.;
    }

    /* Fourier transform */
    kiss_fft(forw,(kiss_fft_cpx *) trace, (kiss_fft_cpx *) fft);
    
    for (i1=0; i1 < nfft; i1++) {
#ifdef SF_HAS_COMPLEX_H
	fft[i1] = cexpf(fft[i1])/nfft;
#else
	fft[i1] = sf_crmul(cexpf(fft[i1]),1./nfft);
#endif
    }

    /* Inverse transform */
    kiss_fft(invs,(const kiss_fft_cpx *) fft, (kiss_fft_cpx *) trace);
}

