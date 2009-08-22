/* Kolmogoroff complex-domain spectral factorization */
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

#include "chelix.h"
/*^*/

#include "kolmog.h"

static kiss_fft_cfg forw, invs;
static int nfft, nk;
static sf_complex *fft1, *fft2;

int xkolmog_init(int n1)
/*< initialize with data length >*/
{
    nfft = n1;
    nk = kiss_fft_next_fast_size(nfft);

    fft1 = sf_complexalloc(nk);   
    fft2 = sf_complexalloc(nk);   

    forw = kiss_fft_alloc(nk,1,NULL,NULL);
    invs = kiss_fft_alloc(nk,0,NULL,NULL);

    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");
    return nk;
}

void xkolmog_close(void)
/*< free allocated storage >*/
{
    free(fft1);
    free(fft2);
    free(forw);
    free(invs);
}

void xkolmog(sf_complex *trace1, sf_complex *trace2)
/*< convert Fourier-domain cross-correlation to minimum-phase >*/ 
{
    int i1;
    const double eps=1.e-32;

    for (i1=0; i1 < nk; i1++) {
#ifdef SF_HAS_COMPLEX_H
	fft1[i1] = clogf(trace1[i1]+eps)/nk;
#else
	fft1[i1] = sf_crmul(clogf(sf_cadd(trace1[i1],sf_cmplx(eps,0.))),
			    1.0/nk);
#endif
    }

    /* Inverse transform */
    kiss_fft(invs,(const kiss_fft_cpx *) fft1, (kiss_fft_cpx *) trace1);

#ifdef SF_HAS_COMPLEX_H
    trace1[0]    *= 0.5; trace2[0]    = trace1[0];
    trace1[nk/2] *= 0.5; trace2[nk/2] = trace1[nk/2];
#else
    trace1[0]    = sf_crmul(trace1[0],   0.5); trace2[0]    = trace1[0];
    trace1[nk/2] = sf_crmul(trace1[nk/2],0.5); trace2[nk/2] = trace1[nk/2];
#endif
    for (i1=1+nk/2; i1 < nk; i1++) {
	trace2[nk-i1] = trace1[i1];
	trace1[i1] = sf_cmplx(0.,0.);
	trace2[i1] = sf_cmplx(0.,0.);
    }

    /* Fourier transform */
    kiss_fft(forw,(const kiss_fft_cpx *) trace1, (kiss_fft_cpx *) fft1);
    kiss_fft(forw,(const kiss_fft_cpx *) trace2, (kiss_fft_cpx *) fft2);
    
    for (i1=0; i1 < nk; i1++) {
#ifdef SF_HAS_COMPLEX_H
	fft1[i1] = cexpf(fft1[i1])/nk;
	fft2[i1] = cexpf(fft2[i1])/nk;
#else
	fft1[i1] = sf_crmul(cexpf(fft1[i1]),1./nk);
	fft2[i1] = sf_crmul(cexpf(fft2[i1]),1./nk);
#endif
    }

    /* Inverse transform */
    kiss_fft(invs,(const kiss_fft_cpx *) fft1, (kiss_fft_cpx *) trace1);
    kiss_fft(invs,(const kiss_fft_cpx *) fft2, (kiss_fft_cpx *) trace2);

    for (i1=0; i1 < nk; i1++) {
	trace2[i1] = conjf(trace2[i1]);
    }
}

void xkolmog_helix(cfilter cross, cfilter fac1, cfilter fac2)
/*< Helix filter factorization >*/
{
    int ik, ih;
    float dw, w;
    
    dw=2.0*SF_PI/nk;

    for (ik=0; ik < nk; ik++) {
	fft1[ik]=sf_cmplx(0.,0.);

	w = dw*ik;

	for (ih=0; ih < cross->nh; ih++) {
#ifdef SF_HAS_COMPLEX_H
	    fft1[ik] += cross->flt[ih]*cexpf(sf_cmplx(0.,cross->lag[ih]*w));
#else
	    fft1[ik] = sf_cadd(fft1[ik],
			       sf_cmul(cross->flt[ih],
				       cexpf(sf_cmplx(0.,cross->lag[ih]*w))));
#endif
	}
#ifdef SF_HAS_COMPLEX_H
	fft1[ik] /= nk;
#else
	fft1[ik] = sf_crmul(fft1[ik],1.0/nk);	
#endif
    }

    xkolmog(fft1,fft2);

    for (ih=0; ih < fac1->nh; ih++) {
	fac1->flt[ih] = fft1[fac1->lag[ih]];
	fac2->flt[ih] = fft2[fac2->lag[ih]];
    }
}
