/* 1-D FFT interface */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

static int n;
static float wt;
static float *ff=NULL;

#ifdef SF_HAS_FFTW
static fftwf_plan cfg=NULL, icfg=NULL;
#else
static kiss_fftr_cfg cfg=NULL, icfg=NULL;
#endif

int fft1_init(int n1  /* input data size */, 
	      int *n2 /* padded data size */)
/*< initialize >*/
{
    int nk;

    nk = kiss_fft_next_fast_size((n1+1)/2)+1;
    n = 2*(nk-1);

    *n2 = n;
    ff = sf_floatalloc(n);

    wt = 1.0/n;

    return nk;
}

void fft1(float *inp      /* [n] */, 
	  sf_complex *out /* [nk] */)
/*< 1-D FFT >*/
{
    int i;

    if (NULL==cfg) {
#ifdef SF_HAS_FFTW
	cfg = fftwf_plan_dft_r2c_1d(n, ff, (fftwf_complex *) out,
				    FFTW_MEASURE);
#else
	cfg  = kiss_fftr_alloc(n,0,NULL,NULL);
#endif
	if (NULL == cfg) sf_error("FFT allocation failure.");
    }

    for (i=0; i < n; i++) {
	ff[i] = inp[i];
    }

#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);
#else
    kiss_fftr (cfg,ff,(kiss_fft_cpx *) out);
#endif
}

void ifft1_allocate(sf_complex *inp /* [nk] */)
/*< allocate inverse transform >*/
{
#ifdef SF_HAS_FFTW
    icfg = fftwf_plan_dft_c2r_1d(n, (fftwf_complex *) inp, ff,
				 FFTW_MEASURE);
    if (NULL == icfg) sf_error("FFT allocation failure.");
#endif
}

void ifft1(float *out      /* [n] */, 
	   sf_complex *inp /* [nk] */)
/*< 1-D inverse FFT >*/
{
    int i;

    if (NULL==icfg) {
#ifdef SF_HAS_FFTW
	sf_error("Call ifft1_allocate first.");
#else
	icfg  = kiss_fftr_alloc(n,1,NULL,NULL);
#endif
	if (NULL == icfg) sf_error("FFT allocation failure.");
    }
    
#ifdef SF_HAS_FFTW
    fftwf_execute(icfg);
#else
    kiss_fftri(icfg,(kiss_fft_cpx *) inp,ff);
#endif
    
    for (i=0; i < n; i++) {
	out[i] = ff[i]*wt;
    }
}

