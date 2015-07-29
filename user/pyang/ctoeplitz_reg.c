/* Regularized Fast Toeplitz matrix inversion: C*mm=dd, C=circulant(c)
   Note: circulant(c) is a kind of teoplitz(c) while satisfying 
   circulant property. 
   C=circulant(c[0],c[1],...,c[n-1])
   Reference: Vogel, Curtis R. Computational methods for inverse problems. 
   Vol. 23. Siam, 2002.	[Chapter 5.2.]
*/

/*
  Copyright (C) 2014 Pengliang Yang, Xi'an Jiaotong University, UT Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

#include "ctoeplitz_reg.h"

#ifdef SF_HAS_FFTW
#include <fftw3.h>


fftwf_plan fft1, ifft1;/* execute plan for FFT and IFFT */
fftwf_complex *tmp;

void ctoeplitz_inv(int n, float mu, sf_complex *c, sf_complex *mm, sf_complex *dd)
/*< linear equation solver with FFT based on Toeplitz structure
  C*mm=dd, C=circulant(c): mm=C^{-1}dd 
  C=F^* diag(c) F 
  ==> C=F^* diag(1./c ) F 
  ==> mm=F^* diag(conj(c)./(c*conj(c) +mu)) F dd
  Here, mu is a stabalizing factor. Setting mu=0 implies no regularization.
  >*/
{
  int i;
  float a;

  tmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n);
  fft1=fftwf_plan_dft_1d(n,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);	
  ifft1=fftwf_plan_dft_1d(n,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
	
  /* dd-->F {dd}*/
  memcpy(tmp, dd, n*sizeof(sf_complex));
  fftwf_execute(fft1);
  memcpy(dd, tmp, n*sizeof(sf_complex));

  /* c-->FFT{c}*/
  memcpy(tmp, c, n*sizeof(sf_complex));
  fftwf_execute(fft1);
  memcpy(c, tmp, n*sizeof(sf_complex));

  /*multiplication in Fourier domain: diag(conj(c)./(c*conj(c)+mu)) F dd*/
  for(i=0; i<n; i++)
    {
      a=c[i]*conjf(c[i]);
      /*dd[i]=(a==0.)?0.:(dd[i]*conjf(c[i])/a);*/
      dd[i]*=conjf(c[i])/(a+mu);
    }

  /* IFFT{FFT{c}.*FFT{dd}/sqrtf(n)}/sqrtf(n)*/
  memcpy(tmp, dd, n*sizeof(sf_complex));
  fftwf_execute(ifft1);
  for(i=0; i<n; i++) mm[i]=tmp[i]/n;
	
}

#endif
