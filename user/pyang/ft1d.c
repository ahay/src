/* 1D fft operator and its adjoint
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
#include <complex.h>
#include <fftw3.h>

#include "ft1d.h"

static int n1;
fftwf_plan fft1, ifft1;/* execute plan for FFT and IFFT */
fftwf_complex *tmp;

void ft1d_init(int n1_)
/*< initialize >*/
{
    	n1=n1_;
    	tmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1);
    	fft1=fftwf_plan_dft_1d(n1,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);	
   	ifft1=fftwf_plan_dft_1d(n1,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
}

void ft1d_lop(bool adj, bool add, int nx, int ny, sf_complex *xx, sf_complex *yy)
/*< linear operator >*/
{
    int i;

    if (n1!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,n1,nx);
    if (n1!=ny) sf_error("%s: size mismatch: %d != %d",__FILE__,n1,ny);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
    for(i=0;i<n1;i++) tmp[i] = adj? yy[i]: xx[i];
  
    if (adj){
 	fftwf_execute(fft1);	
    	for(i=0;i<n1;i++) xx[i]+=tmp[i]/sqrtf(n1);
    }else{
	fftwf_execute(ifft1);
    	for(i=0;i<n1;i++) yy[i]+=tmp[i]/sqrtf(n1);
    }
}

void ft1d_close(void)
/*< free allocated variables>*/
{
    fftwf_free(tmp);
    fftwf_destroy_plan(fft1);
    fftwf_destroy_plan(ifft1);
}

