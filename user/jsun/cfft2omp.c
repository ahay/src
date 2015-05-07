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

static int n1, n2, nk;
static float wt;
static sf_complex **cc, **trace2;

static kiss_fft_cfg *cfg1, *icfg1, *cfg2, *icfg2;
static kiss_fft_cpx **tmp, **ctrace2;


int cfft2_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny   /* input data size */, 
	       int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{

  int i2,i;
  int nth=1;

  nk = n1 = kiss_fft_next_fast_size(nx*pad1);
  n2 = kiss_fft_next_fast_size(ny);

  cc = sf_complexalloc2(n1,n2);
    
#ifdef _OPENMP
#pragma omp parallel
  {nth = omp_get_num_threads();}
#endif

  cfg1  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
  icfg1 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
  cfg2  = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));
  icfg2 = (kiss_fft_cfg *) sf_alloc(nth,sizeof(kiss_fft_cfg));

  for (i=0; i < nth; i++) {
    cfg1[i] = kiss_fft_alloc(n1,0,NULL,NULL);
    icfg1[i]= kiss_fft_alloc(n1,1,NULL,NULL);
    cfg2[i] = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2[i]= kiss_fft_alloc(n2,1,NULL,NULL);
  }

  trace2 = sf_complexalloc2(n2,nth);
  ctrace2= (kiss_fft_cpx **) trace2;

  tmp =    (kiss_fft_cpx **) sf_alloc(n2,sizeof(*tmp));
  tmp[0] = (kiss_fft_cpx *)  sf_alloc(nk*n2,sizeof(kiss_fft_cpx));
#ifdef _OPENMP
#pragma omp parallel for private(i2) default(shared)
#endif
  for (i2=0; i2 < n2; i2++) {
    tmp[i2] = tmp[0]+i2*nk;
  }

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
  int ith=0;

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

#ifdef _OPENMP
#pragma omp parallel for private(i2,ith) default(shared)
#endif
  for (i2=0; i2 < n2; i2++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(cfg1[ith],(kiss_fft_cpx *) cc[i2],tmp[i2],1);
  }

#ifdef _OPENMP
#pragma omp parallel for private(i1,i2,ith) default(shared)
#endif
  for (i1=0; i1 < nk; i1++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(cfg2[ith],tmp[0]+i1,ctrace2[ith],nk);
    for (i2=0; i2<n2; i2++) {
      out[i2*nk+i1] = trace2[ith][i2];
    }
  }

}

void icfft2(sf_complex *out /* [n1*n2] */, 
	    sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
  int i1, i2;
  int ith=0;

  //  sf_warning("OK-fan");

#ifdef _OPENMP
#pragma omp parallel for private(i1,i2,ith) default(shared)
#endif
  for (i1=0; i1 < nk; i1++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(icfg2[ith],(kiss_fft_cpx *) (inp+i1),ctrace2[ith],nk);
    for (i2=0; i2<n2; i2++) {
      tmp[i2][i1] = ctrace2[ith][i2];
    }
  }

#ifdef _OPENMP
#pragma omp parallel for private(i2,ith) default(shared)
#endif
  for (i2=0; i2 < n2; i2++) {
#ifdef _OPENMP
    ith = omp_get_thread_num();
#endif
    kiss_fft_stride(icfg1[ith],tmp[i2],(kiss_fft_cpx *) cc[i2],1);
  }
    
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
/*< clean up fft >*/
{
  /* make sure everything is back to its pristine state */
  int nth=1,i;

#ifdef _OPENMP
#pragma omp parallel
  {nth = omp_get_num_threads();}
#endif
  for (i=0; i < nth; i++) {
    free(cfg1[i]);   cfg1[i]=NULL;
    free(icfg1[i]); icfg1[i]=NULL;
    free(cfg2[i]);   cfg2[i]=NULL;
    free(icfg2[i]); icfg2[i]=NULL;
  }

  free(*cc); free(cc);
  free(*tmp); free(tmp);
  free(*ctrace2); free(ctrace2);
}
