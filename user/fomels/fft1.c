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

#ifdef SF_HAS_FFTW
static fftwf_plan cfg=NULL, *icfg=NULL;
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
    return nk;
}

void ifft1_allocate(int ni /* number of inverse transforms */)
/*< allocate inverse transforms >*/
{
#ifdef SF_HAS_FFTW
    int i;

    icfg = (fftwf_plan*) sf_alloc(ni,sizeof(fftwf_plan));
    for (i=0; i < ni; i++) {
	icfg[i] = NULL;
    }
#endif
}

void fft1(float *inp      /* [n] */, 
	  sf_complex *out /* [nk] */)
/*< 1-D FFT >*/
{
    if (NULL==cfg) {
#ifdef SF_HAS_FFTW
	cfg = fftwf_plan_dft_r2c_1d(n,
				    inp, (fftwf_complex *) out,
				    FFTW_ESTIMATE);
#else
	cfg  = kiss_fftr_alloc(n,0,NULL,NULL);
#endif
	if (NULL == cfg) sf_error("FFT allocation failure.");
    }
#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);
#else
    kiss_fftr (cfg,inp,(kiss_fft_cpx *) out);
#endif
}


void ifft1(int i           /* transform number */,
	   float *out      /* [n] */, 
	   sf_complex *inp /* [nk] */)
/*< 1-D inverse FFT >*/
{
#ifdef SF_HAS_FFTW
    if (NULL==icfg[i]) {
	icfg[i] = fftwf_plan_dft_c2r_1d(n, (fftwf_complex *) inp, out,
					FFTW_ESTIMATE);
	if (NULL == icfg[i]) sf_error("FFT allocation failure.");
    }
    fftwf_execute(icfg[i]);
#else
    if (NULL==icfg) {
	icfg  = kiss_fftr_alloc(n,1,NULL,NULL);
	if (NULL == icfg) sf_error("FFT allocation failure.");
    }
    kiss_fftri(icfg,(kiss_fft_cpx *) inp,out);
#endif
}

