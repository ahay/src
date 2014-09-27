/* discrete linear chirp transfrom (DLCT)
   To make the adjoint as same as the inverse, I normalized the forward 
   transform of DLCT with a factor sqrt(N*L).
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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

  Reference: Alkishriwo, Osama A., and Luis F. Chaparro. "A discrete 
  linear chirp transform (DLCT) for data compression." Information 
  Science, Signal Processing and their Applications (ISSPA), 2012 
  11th International Conference on. IEEE, 2012.

  We can think of DLCT as a complex-valued Gabor transform! Complex-valued 
  version does not introduce amplitude scaling!
*/

#include <rsf.h>
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "dlct.h"

#ifdef SF_HAS_FFTW
#include <fftw3.h>

static int N, L;
static float C;
fftwf_complex *tmp;
fftwf_plan fft, ifft;

void dlct_init( int N_ 	/* length of the signal */,
		int L_	/* length of freq-instaneous freq */, 
		float C_/* step size for freq-inst freq */)
/*< initialize DLCT transform >*/
{
  N=N_;
  L=L_;
  C=C_;
  tmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);
  fft=fftwf_plan_dft_1d(N, tmp, tmp,FFTW_FORWARD,FFTW_MEASURE);
  ifft=fftwf_plan_dft_1d(N, tmp, tmp,FFTW_BACKWARD,FFTW_MEASURE);
}


void dlct_close(void)
/*< free allocated storage >*/
{
  fftwf_destroy_plan(fft);
  fftwf_destroy_plan(ifft);
  fftwf_free(tmp);
}

void dlct_lop(bool adj, bool add, int nm, int nd, sf_complex *mm, sf_complex *dd)
/*< DLCT linear operator >*/
{

  int n,k,l;
  if(nm!=N*L || nd!=N) sf_error("datasize mismatch!");
  sf_cadjnull (adj, add,  nm, nd, mm, dd); 

  if(adj){// forward transform: dd-->mm (data-->model)
    for(l=-L/2;l<L/2;l++){
#ifdef _OPENMP
#pragma omp parallel for default(none) private(n) shared(tmp,dd,C,N,l)
#endif	
      for(n=0;n<N;n++) tmp[n]=dd[n]*cexpf(-I*2.*SF_PI*C*l*n*n/N);
      fftwf_execute(fft);
#ifdef _OPENMP
#pragma omp parallel for default(none) private(k) shared(mm,l,L,N,tmp)
#endif	
      for(k=0;k<N;k++) mm[k+(l+L/2)*N]+=tmp[k]/sqrtf(N*L);
    }
  }else{// inverse transform: mm-->dd (model-->data)
    for(l=-L/2;l<L/2;l++){
#ifdef _OPENMP
#pragma omp parallel for default(none) private(k) shared(tmp,mm,N,L,l)
#endif	
      for(k=0;k<N;k++) tmp[k]=mm[k+(l+L/2)*N];
      fftwf_execute(ifft);
#ifdef _OPENMP
#pragma omp parallel for default(none) private(n) shared(dd,tmp,C,l,N,L)
#endif	
      for(n=0;n<N;n++) dd[n]+=tmp[n]*cexpf(I*2.*SF_PI*C*l*n*n/N)/sqrtf(N*L);
    }
  }
}

#endif
