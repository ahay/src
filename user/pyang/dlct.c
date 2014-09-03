/* discrete linear chirp transfrom (DLCT)
   Note: In my implementation:, to make the adjoint as same as the inverse,
   I normalized the forward transform of DLCT with a factor sqrt(N*L).
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "dlct.h"

#ifdef SF_HAS_FFTW
#include <fftw3.h>

void forward_dlct(int N 	/* length of the signal */,
		  int L		/* length of freq-instaneous freq */, 
		  float C		/* step size for freq-inst freq */,
		  float *d	/* input [N] signal float or complex */,
		  sf_complex *Sc	/* output[N*L] DLCT coefficients */ )
/*< forward DLCT >*/
{
  int n,k,l;
  fftwf_complex *p,*q;
  fftwf_plan fft;

  p=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);
  q=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);

  fft=fftwf_plan_dft_1d(N,p, q,FFTW_FORWARD,FFTW_MEASURE);

  for(l=-L/2;l<L/2;l++){

#ifdef _OPENMP
#pragma omp parallel for default(none) private(n) shared(p,d,C,N,l)
#endif	
    for(n=0;n<N;n++){
#ifdef SF_HAS_COMPLEX_H
      p[n]=d[n]*cexpf(-I*2*SF_PI*C*l*n*n/N);
#else
      p[n]=sf_crmul(d[n],cexpf(sf_cmplx(0, -2*SF_PI*C*l*n*n/N)));
#endif
    }	
    fftwf_execute(fft);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(k) shared(Sc,l,L,N,q)
#endif	
    for(k=0;k<N;k++)	Sc[k+(l+L/2)*N]=q[k]/sqrtf(N*L);
  }
  fftwf_destroy_plan(fft);
  fftwf_free(p);
  fftwf_free(q);
}


void inverse_dlct(int N 	/* length of the signal */,
		  int L		/* length of freq-instaneous freq */, 
		  float C		/* step size for freq-inst freq */,
		  float *d	/* output [N] signal,float or complex */,
		  sf_complex *Sc	/* input[N*L] DLCT coefficients */ )
/*< inverse DLCT >*/
{  
  int l,k,n;
  fftwf_complex *p,*q;
  fftwf_plan ifft;
  p=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);
  q=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);

  ifft=fftwf_plan_dft_1d(N,p, q,FFTW_BACKWARD,FFTW_MEASURE);

  memset(d,0,N*sizeof(float));
  for(l=-L/2;l<L/2;l++){

#ifdef _OPENMP
#pragma omp parallel for default(none) private(k) shared(p,Sc,N,L,l)
#endif	
    for(k=0;k<N;k++) 	p[k]=Sc[k+(l+L/2)*N];

    fftwf_execute(ifft);

#ifdef _OPENMP
#pragma omp parallel for default(none) private(n) shared(d,q,C,l,N,L)
#endif	
    for(n=0;n<N;n++){	
#ifdef SF_HAS_COMPLEX_H	
      d[n]+=crealf(q[n]*cexpf(I*2*SF_PI*C*l*n*n/N)/sqrtf(N*L));
#else
      d[n]=sf_cadd(d[n],crealf(q[n]*cexpf(sf_cmplx(0, 2*SF_PI*C*l*n*n/N))/sqrtf(N*L)));
#endif
    }
  }

  fftwf_destroy_plan(ifft);
  fftwf_free(p);
  fftwf_free(q);
}

#endif
