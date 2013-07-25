/* Complex 2-D FFT interface, complex to complex fft and ifft */
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

#include <fftw3.h>


static int n1, n2, nk;
static float wt;

static sf_complex **cc;


static fftwf_plan cfg=NULL, icfg=NULL;

int wcfft2_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny   /* input data size */, 
	       int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{
#ifdef _OPENMP
    fftw_init_threads();
    //sf_warning("Using threaded FFTW3! \n");
#endif

    nk = n1 = kiss_fft_next_fast_size(nx*pad1);
    
    n2 = kiss_fft_next_fast_size(ny);

    cc = sf_complexalloc2(n1,n2);
    
    *nx2 = n1;
    *ny2 = n2;
	
    wt =  1.0/(n1*n2);
	
    return (nk*n2);
}

void wcfft2(sf_complex *inp /* [n1*n2] */, 
	   sf_complex *out /* [nk*n2] */)
/*< 2-D FFT >*/
{
    int i1, i2;

#ifdef _OPENMP
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    if (NULL==cfg) {
      cfg = fftwf_plan_dft_2d(n2,n1,
			      (fftwf_complex *) cc[0], 
			      (fftwf_complex *) out,
			      FFTW_FORWARD, FFTW_MEASURE);
      if (NULL == cfg) sf_error("FFTW failure.");
    }

    /* FFT centering */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
		cc[i2][i1] = ((i2%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1];
	}
    }

    fftwf_execute(cfg);

}

void iwcfft2_allocate(sf_complex *inp /* [nk*n2] */)
/*< allocate inverse transform >*/
{
    icfg = fftwf_plan_dft_2d(n2,n1,
			     (fftwf_complex *) inp, 
			     (fftwf_complex *) cc[0],
			     FFTW_BACKWARD, FFTW_MEASURE);
    if (NULL == icfg) sf_error("FFTW failure.");
}

void iwcfft2(sf_complex *out /* [n1*n2] */, 
	    sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
    int i1, i2;

    fftwf_execute(icfg);

    /* FFT centering and normalization*/
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

void wcfft2_finalize()
/*< clean up fftw >*/
{
#ifdef _OPENMP
    fftw_cleanup_threads();
#endif
}
