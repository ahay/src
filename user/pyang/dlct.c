/* discrete linear chirp transfrom (DLCT)
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
*/

#include <rsf.h>
#include <complex.h>
#include <fftw3.h>

#include "dlct.h"

#ifdef _OPENMP
#include <omp.h>
#endif


void forward_dlct(int N 	/* length of the signal */,
		int L		/* length of freq-instaneous freq */, 
		float C		/* step size for freq-inst freq */,
		float *d		/* input [N] signal float or complex */,
		sf_complex *Sc	/* output[N*L] DLCT coefficients */ )
/*< forward DLCT >*/
{
    fftwf_complex *p,*q;
    fftwf_plan fft;

    p=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);
    q=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);

    fft=fftwf_plan_dft_1d(N,p, q,FFTW_FORWARD,FFTW_MEASURE);

    for(int l=-L/2;l<L/2;l++){
	for(int n=0;n<N;n++){
	  p[n]=d[n]*cexpf(sf_cmplx(0, -2*SF_PI*C*l*n*n/N));
	}	
	fftwf_execute(fft);

	for(int k=0;k<N;k++)	Sc[k+(l+L/2)*N]=q[k];
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
    fftwf_complex *p,*q;
    fftwf_plan ifft;
    p=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);
    q=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*N);

    ifft=fftwf_plan_dft_1d(N,p, q,FFTW_BACKWARD,FFTW_MEASURE);

    for(int n=0;n<N;n++)  {
	d[n]=0.0;
    }
    for(int l=-L/2;l<L/2;l++){
	for(int k=0;k<N;k++) 	p[k]=Sc[k+(l+L/2)*N];

	fftwf_execute(ifft);
	for(int n=0;n<N;n++){		
	  d[n]+=crealf(q[n]*cexpf(sf_cmplx(0, 2*SF_PI*C*l*n*n/N))/(N*L));
	}
    }

    fftwf_destroy_plan(ifft);
    fftwf_free(p);
    fftwf_free(q);
}
