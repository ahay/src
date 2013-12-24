/* 3D mask_fft linear operator interface */

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

#include "maskft3d.h"

static int n1, n2, n3;
float *mask;
fftwf_plan fft3, ifft3;/* execute plan for FFT and IFFT */
sf_complex *tmp;

void maskft3d_init(int n1_, int n2_, int n3_, float *mask_)
/*< initialize >*/
{
    	n1=n1_;
    	n2=n2_;
	n3=n3_;
	mask=(float*)malloc(n2*n3);
	for(int i=0;i<n2*n3;i++) mask[i]=mask_[i];
    	tmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    	fft3=fftwf_plan_dft_3d(n1,n2,n3,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);	
   	ifft3=fftwf_plan_dft_3d(n1,n2,n3,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
}

void maskft3d_lop (bool adj, bool add, int nx, int ny, sf_complex *xx, sf_complex *yy)
/*< linear operator >*/
{
    if (n1*n2*n3!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__,n1*n2*n3,nx);
    if (n1*n2*n3!=ny) sf_error("%s: size mismatch: %d != %d",__FILE__,n1*n2*n3,ny);

    sf_cadjnull (adj, add, nx, ny, xx, yy);// if adj, null(xx); else null(yy)
    for(int i=0;i<n1*n2*n3;i++) tmp[i] = adj? yy[i]: xx[i];
  
    if (adj){
 	fftwf_execute(ifft3);	

	for(int i3=0;i3<n3;i3++)	
	for(int i2=0;i2<n2;i2++)
	{
		if (mask[i2+i3*n2]){			
		    for(int i1=0; i1<n1; i1++) xx[i1]+=tmp[i1];
		}
	}
    }else{
	fftwf_execute(fft3);

	for(int i3=0;i3<n3;i3++)	
	for(int i2=0;i2<n2;i2++)
	{
		if (mask[i2+i3*n2]){			
		    for(int i1=0; i1<n1; i1++) yy[i1]+=tmp[i1];
		}
	}
    }
}

void maskft3d_close(void)
/*< free allocated variables >*/
{
    free(mask);
    fftwf_free(tmp);
    fftwf_destroy_plan(fft3);
    fftwf_destroy_plan(ifft3);
}

