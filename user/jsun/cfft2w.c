/* Complex 2-D FFT interface, complex to complex fft and ifft with threaded FFTW3, designed for NSPS operator */
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
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

static int n1, n2, nk;
static float wt;
static sf_complex **cc;
#ifdef SF_HAS_FFTW
static sf_complex **dd;
static fftwf_plan cfg=NULL, icfg=NULL;
#else
static kiss_fft_cfg cfg1, icfg1, cfg2, icfg2;
static kiss_fft_cpx **tmp, *ctrace2;
static sf_complex *trace2;
#endif

int cfft2_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny   /* input data size */, 
	       int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{

#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_init_threads();
    if (false)
      sf_warning("Using threaded FFTW3! \n");
    fftwf_plan_with_nthreads(omp_get_max_threads());
#else
    if (false)
      sf_warning("Using FFTW3! \n");
#endif
#else
    if (false)
      sf_warning("Using KissFFT! \n");
#endif

#ifndef SF_HAS_FFTW
    int i2;
#endif

    nk = n1 = kiss_fft_next_fast_size(nx*pad1);
    n2 = kiss_fft_next_fast_size(ny);

    cc = sf_complexalloc2(n1,n2);
    
#ifdef SF_HAS_FFTW
    dd = sf_complexalloc2(nk,n2);

    cfg = fftwf_plan_dft_2d(n2,n1,
			    (fftwf_complex *) cc[0], 
			    (fftwf_complex *) dd[0],
			    FFTW_FORWARD, FFTW_MEASURE);

    icfg = fftwf_plan_dft_2d(n2,n1,
			     (fftwf_complex *) dd[0], 
			     (fftwf_complex *) cc[0],
			     FFTW_BACKWARD, FFTW_MEASURE);

    if (NULL == cfg || NULL == icfg) sf_error("FFTW failure.");
#else
    cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
    icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);

    cfg2  = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(n2,1,NULL,NULL);
 	
    /* temporary array */

    tmp =    (kiss_fft_cpx **) sf_alloc(n2,sizeof(*tmp));
    tmp[0] = (kiss_fft_cpx *)  sf_alloc(nk*n2,sizeof(kiss_fft_cpx));
#ifdef _OPENMP
#pragma omp parallel for private(i2) default(shared)
#endif
    for (i2=0; i2 < n2; i2++) {
	tmp[i2] = tmp[0]+i2*nk;
    }
	
    trace2 = sf_complexalloc(n2);
    ctrace2 = (kiss_fft_cpx *) trace2;
#endif

    *nx2 = n1;
    *ny2 = n2;
	
    wt =  1.0/(n1*n2);
	
    return (nk*n2);
}

void cfft2(sf_complex *inp /* [n1*n2] */, 
	   sf_complex *out /* [nk*n2] */)
/*< 2-D FFT >*/
{
    int i1, i2;

    /* FFT centering */
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		cc[i2][i1] = ((i2%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1];
#else
		cc[i2][i1] = ((i2%2==0)==(i1%2==0))? inp[i2*n1+i1]:sf_cneg(inp[i2*n1+i1]);
#endif
	}
    }

#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<nk; i1++) {
	    out[i2*nk+i1]=dd[i2][i1];
	}
    }
#else
    for (i2=0; i2 < n2; i2++) {
	kiss_fft_stride(cfg1,(kiss_fft_cpx *) cc[i2],tmp[i2],1);
    }
    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(cfg2,tmp[0]+i1,ctrace2,nk);
	for (i2=0; i2<n2; i2++) {
	    out[i2*nk+i1] = trace2[i2];
	}
    }
#endif
}

void icfft2_allocate(sf_complex *inp /* [nk*n2] */)
/*< allocate inverse transform >*/
{
  /*kept for backward compatibility*/
}

void icfft2(sf_complex *out /* [n1*n2] */, 
	    sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<nk; i1++) {
	    dd[i2][i1]=inp[i2*nk+i1];
	}
    }
    fftwf_execute(icfg);
#else
    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(icfg2,(kiss_fft_cpx *) (inp+i1),ctrace2,nk);
	for (i2=0; i2<n2; i2++) {
	    tmp[i2][i1] = ctrace2[i2];
	}
    }
    for (i2=0; i2 < n2; i2++) {
	kiss_fft_stride(icfg1,tmp[i2],(kiss_fft_cpx *) cc[i2],1);
    }
#endif
    
    /* FFT centering and normalization*/
#ifdef _OPENMP
#pragma omp parallel for private(i2,i1) default(shared)
#endif
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		out[i2*n1+i1] = (((i2%2==0)==(i1%2==0))? wt:-wt) * cc[i2][i1];
#else
		out[i2*n1+i1] = sf_crmul(cc[i2][i1],(((i2%2==0)==(i1%2==0))? wt:-wt));
#endif
	}
    }
}

void cfft2_finalize()
/*< clean up fftw >*/
{
/* make sure everything is back to its pristine state */
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_cleanup_threads();
#endif
    fftwf_destroy_plan(cfg);
    fftwf_destroy_plan(icfg);
    fftwf_cleanup();
    cfg=NULL;
    icfg=NULL;
    free(*dd);
    free(dd);
#else
    free(cfg1); cfg1=NULL;
    free(icfg1); icfg1=NULL;
    free(cfg2); cfg2=NULL;
    free(icfg2); icfg2=NULL;
#endif
    free(*cc);
    free(cc);
}
