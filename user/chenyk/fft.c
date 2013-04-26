/* 1D fourier tranform unitlity */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include<rsf.h>
#include"fft.h"

void fft(float *xx /*input data*/, 
		    kiss_fft_cpx *xxf  /*spectrum*/, 
	 	    int n1 /*length of array*/)
/*< forward fourier tranform  >*/
{

   int nfft;
   kiss_fftr_cfg forw;
   /* determine frequency sampling (for real to complex FFT) */
   nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
   forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
   if(NULL==forw )
    	sf_error("%s: KISS FFT allocation problem", __FILE__);
   kiss_fftr(forw, xx, xxf);
}

void ifft(kiss_fft_cpx *xxf /*spectrum*/, 
		    float *xx  /*output data*/, 
	 	    int n1 /*length of array*/)
/*< inverse fourier tranform >*/
{

   int nfft,nw,i;
   kiss_fftr_cfg invs;
   /* determine frequency sampling (for real to complex FFT) */
   nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
   nw=nfft/2+1;
   invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
   if(NULL==invs)
   	sf_error("%s: KISS FFT allocation problem", __FILE__);
   for(i=0;i<nw;i++)
        xxf[i]=sf_crmul(xxf[i],1.0/nfft); 
   kiss_fftri(invs,xxf,xx);
}
