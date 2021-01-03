/* 2-D FFT interface with threaded FFTW3 */
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

static bool cmplx;
static int n1, n2, nk;
static float wt;

static float **ff=NULL;
static sf_complex **cc=NULL,*dd=NULL;

#ifdef SF_HAS_FFTW
static fftwf_plan cfg=NULL, icfg=NULL;
#else
static kiss_fftr_cfg cfg=NULL, icfg=NULL;
static kiss_fft_cfg cfg1=NULL, icfg1=NULL, cfg2=NULL, icfg2=NULL;
static kiss_fft_cpx **tmp=NULL, *ctrace2=NULL;
static sf_complex *trace2=NULL;
#endif

int fft2_init(bool cmplx1        /* if complex transform */,
	      int pad1           /* padding on the first axis */,
	      int nx,   int ny   /* input data size */, 
	      int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{
#ifdef SF_HAS_FFTW
#ifdef _OPENMP
    fftwf_init_threads();
    sf_warning("Using threaded FFTW3!\n");
    fftwf_plan_with_nthreads(omp_get_max_threads());
#endif
#else
   int i2;
#endif
	
    cmplx = cmplx1;
	
    if (cmplx) {
	nk = n1 = kiss_fft_next_fast_size(nx*pad1);
		
#ifndef SF_HAS_FFTW
	cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
	icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);
#endif
    } else {
	nk = kiss_fft_next_fast_size(pad1*(nx+1)/2)+1;
	n1 = 2*(nk-1);
		
#ifndef SF_HAS_FFTW
	cfg  = kiss_fftr_alloc(n1,0,NULL,NULL);
	icfg = kiss_fftr_alloc(n1,1,NULL,NULL);
#endif
    }
		
    n2 = kiss_fft_next_fast_size(ny);

    if (cmplx) {
	cc = sf_complexalloc2(n1,n2);
    } else {
	ff = sf_floatalloc2(n1,n2);
    }
    dd = sf_complexalloc(nk*n2);
	
#ifndef SF_HAS_FFTW
    cfg2  = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(n2,1,NULL,NULL);
 	
    tmp =    (kiss_fft_cpx **) sf_alloc(n2,sizeof(*tmp));
    tmp[0] = (kiss_fft_cpx *)  sf_alloc(nk*n2,sizeof(kiss_fft_cpx));
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

void fft2(float *inp      /* [n1*n2] */, 
	  sf_complex *out /* [nk*n2] */)
/*< 2-D FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    if (NULL==cfg) {
	cfg = cmplx? 
	    fftwf_plan_dft_2d(n2,n1,
			      (fftwf_complex *) cc[0],
			      (fftwf_complex *) dd,
			      FFTW_FORWARD, FFTW_MEASURE):
	    fftwf_plan_dft_r2c_2d(n2,n1,
				  ff[0], (fftwf_complex *) dd,
				  FFTW_MEASURE);
	if (NULL == cfg) sf_error("FFTW failure.");
    }
#endif

    /* FFT centering */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    if (cmplx) {
		cc[i2][i1] = sf_cmplx(((i2%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1],0.);
	    } else {
		ff[i2][i1] = (i2%2)? -inp[i2*n1+i1]:inp[i2*n1+i1];
	    }
	}
    }
    
#ifdef SF_HAS_FFTW
    fftwf_execute(cfg);

    for (i1=0; i1 < nk*n2; i1++)
      out[i1] = dd[i1];
#else	
    for (i2=0; i2 < n2; i2++) {
	if (cmplx) {
	    kiss_fft_stride(cfg1,(kiss_fft_cpx *) cc[i2],tmp[i2],1);
	} else {
	    kiss_fftr (cfg,ff[i2],tmp[i2]);
	}
    }
	
    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(cfg2,tmp[0]+i1,ctrace2,nk);
	for (i2=0; i2<n2; i2++) {
	    out[i2*nk+i1] = trace2[i2];
	}
    }
#endif
}

void ifft2_allocate(sf_complex *inp /* [nk*n2] */)
/*< allocate inverse transform >*/
{
  /* for backward compatibility */
}

void ifft2(float *out      /* [n1*n2] */, 
	   sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
    int i1, i2;

#ifdef SF_HAS_FFTW
    if (NULL==icfg) {
      icfg = cmplx? 
	fftwf_plan_dft_2d(n2,n1,
			  (fftwf_complex *) dd, 
			  (fftwf_complex *) cc[0],
			  FFTW_BACKWARD, FFTW_MEASURE):
	fftwf_plan_dft_c2r_2d(n2,n1,
			      (fftwf_complex *) dd, ff[0],
			      FFTW_MEASURE);
      if (NULL == icfg) sf_error("FFTW failure.");
    }
#endif

#ifdef SF_HAS_FFTW
    for (i1=0; i1 < nk*n2; i1++)
      dd[i1] = inp[i1];

    fftwf_execute(icfg);
#else
    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(icfg2,(kiss_fft_cpx *) (inp+i1),ctrace2,nk);
		
	for (i2=0; i2<n2; i2++) {
	    tmp[i2][i1] = ctrace2[i2];
	}
    }
    for (i2=0; i2 < n2; i2++) {
	if (cmplx) {
	    kiss_fft_stride(icfg1,tmp[i2],(kiss_fft_cpx *) cc[i2],1);
	} else {
	    kiss_fftri(icfg,tmp[i2],ff[i2]);
	}
    }
#endif
    
    /* FFT centering and normalization */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    if (cmplx) {
		out[i2*n1+i1] = (((i2%2==0)==(i1%2==0))? wt:-wt) * crealf(cc[i2][i1]);
	    } else {
		out[i2*n1+i1] = (i2%2? -wt: wt)*ff[i2][i1];
	    }
	}
    }
}

void fft2_finalize()
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
#else
    if (NULL != cfg) { free(cfg); cfg=NULL; }
    if (NULL != icfg) { free(icfg); icfg=NULL; }
    if (NULL != cfg1) { free(cfg1); cfg1=NULL; }
    if (NULL != icfg1) { free(icfg1); icfg1=NULL; }
    if (NULL != cfg2) { free(cfg2); cfg2=NULL; }
    if (NULL != icfg2) { free(icfg2); icfg2=NULL; }
    if (NULL != tmp) { free(*tmp); free(tmp); tmp=NULL; }
    if (NULL != trace2) { free(trace2); trace2=NULL; }
#endif
    if (cmplx) {
      if (NULL != cc) { free(*cc); free(cc); cc=NULL; }
    } else {
      if (NULL != ff) { free(*ff); free(ff); ff=NULL; }
    }
    if (NULL != dd) { free(dd); dd=NULL; }
}

