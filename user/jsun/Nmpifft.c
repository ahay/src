/* Complex 2-D wave propagation (with multi-threaded FFTW3)*/
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <omp.h>
#include <fftw3-mpi.h>

static int n1, n2, nk;
static float wt;
static sf_complex *cc,*dd;
static fftwf_plan cfg=NULL, icfg=NULL;
static ptrdiff_t alloc_local, local_n0, local_0_start;

int threads_ok;

int cfft2_init(int pad1           /* padding on the first axis */,
	       int nx,   int ny   /* input data size */, 
	       int *nx2, int *ny2 /* padded data size */,
               int *n_local, int *o_local /* local size & start */,
               MPI_Comm comm)
/*< initialize >*/
{
  if (threads_ok) threads_ok = fftwf_init_threads();

  fftwf_mpi_init();

  if (false)
    sf_warning("Using threaded FFTW3! \n");
  if (threads_ok)
    fftwf_plan_with_nthreads(omp_get_max_threads());

  nk = n1 = kiss_fft_next_fast_size(nx*pad1);
  n2 = kiss_fft_next_fast_size(ny);

  alloc_local = fftwf_mpi_local_size_2d(n2, n1, comm, &local_n0, &local_0_start);

  //cc = sf_complexalloc2(n1,n2);
  //dd = sf_complexalloc2(nk,n2);
  cc = sf_complexalloc(alloc_local);
  dd = sf_complexalloc(alloc_local);

  cfg = fftwf_mpi_plan_dft_2d(n2,n1,
                              (fftwf_complex *) cc,
                              (fftwf_complex *) dd,
                              comm,
                              FFTW_FORWARD, FFTW_MEASURE);

  icfg = fftwf_mpi_plan_dft_2d(n2,n1,
                               (fftwf_complex *) dd, 
                               (fftwf_complex *) cc,
                               comm,
                               FFTW_BACKWARD, FFTW_MEASURE);

  if (NULL == cfg || NULL == icfg) sf_error("FFTW failure.");

  *nx2 = n1;
  *ny2 = n2;
  *n_local = (int) local_n0;
  *o_local = (int) local_0_start;
	
  wt =  1.0/(n1*n2);
	
  return (nk*n2);
}

void cfft2(sf_complex *inp /* [n1*n2] */, 
	   sf_complex *out /* [nk*n2] */)
/*< 2-D FFT >*/
{
  int i1, i2;

  /* FFT centering */
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
      cc[i2*n1+i1] = (((i2+local_0_start)%2==0)==(i1%2==0))? inp[i2*n1+i1]:-inp[i2*n1+i1];
#else
      cc[i2*n1+i1] = (((i2+local_0_start)%2==0)==(i1%2==0))? inp[i2*n1+i1]:sf_cneg(inp[i2*n1+i1]);
#endif
    }
  }

  fftwf_execute(cfg);
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<nk; i1++) {
      out[i2*nk+i1]=dd[i2*nk+i1];
    }
  }
}

void icfft2(sf_complex *out /* [n1*n2] */, 
	    sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
  int i1, i2;

#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<nk; i1++) {
      dd[i2*nk+i1]=inp[i2*nk+i1];
    }
  }
  fftwf_execute(icfg);
    
  /* FFT centering and normalization*/
#pragma omp parallel for private(i2,i1) default(shared)
  for (i2=0; i2<local_n0; i2++) {
    for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
      out[i2*n1+i1] = ((((i2+local_0_start)%2==0)==(i1%2==0))? wt:-wt) * cc[i2*n1+i1];
#else
      out[i2*n1+i1] = sf_crmul(cc[i2*n1+i1],((((i2+local_0_start)%2==0)==(i1%2==0))? wt:-wt));
#endif
    }
  }
}

void cfft2_finalize()
/*< clean up fftw >*/
{
  /* make sure everything is back to its pristine state */

  fftwf_destroy_plan(cfg);
  fftwf_destroy_plan(icfg);
  fftwf_mpi_cleanup();
  //fftwf_cleanup_threads();
  //fftwf_cleanup();
  cfg=NULL;
  icfg=NULL;
  free(dd);
  free(cc);
}

