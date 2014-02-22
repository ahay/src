/* 3D fft operator */

/*
  Copyright (C) 2013 Pengliang Yang, Xi'an Jiaotong University, UT Austin
  
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

#include "ft3d.h"

static int n1, n2, n3;
fftwf_plan fft3, ifft3;/* execute plan for FFT and IFFT */
sf_complex *tmp;

void ft3d_init(int n1_, int n2_, int n3_)
/*< initialize >*/
{
    	n1=n1_;
    	n2=n2_;
	n3=n3_;
    	tmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    	fft3=fftwf_plan_dft_3d(n1,n2,n3,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);	
   	ifft3=fftwf_plan_dft_3d(n1,n2,n3,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
}

void ft3d_lop (bool adj, bool add, int nx, int ny, sf_complex *xx, sf_complex *yy)
/*< linear operator >*/
{
    int i;

    if (n1*n2*n3!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,n1*n2*n3,nx);
    if (n1*n2*n3!=ny) sf_error("%s: size mismatch: %d != %d",__FILE__,n1*n2*n3,ny);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
    for(i=0;i<n1*n2*n3;i++) tmp[i] = adj? yy[i]: xx[i];
  
    if (adj){
 	fftwf_execute(fft3);	
    	for(i=0;i<n1*n2*n3;i++) xx[i]+=tmp[i]/sqrtf(n1*n2*n3);
    }else{
	fftwf_execute(ifft3);
    	for(i=0;i<n1*n2*n3;i++) yy[i]+=tmp[i]/sqrtf(n1*n2*n3);
    }
}

void ft3d_close(void)
/*< free allocated variables>*/
{
    fftwf_free(tmp);
    fftwf_destroy_plan(fft3);
    fftwf_destroy_plan(ifft3);
}

