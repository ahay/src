/* Frequency domain filter */
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

#include <rsf.h>
#include "ffilter.h"

static int iw, nfft, nw;
static kiss_fft_cpx *cdata;
static float *shape, dw;
static kiss_fftr_cfg forw, invs;

void trapezoid_init(int n /*length of time series*/, 
		    float fn /*nyquist frequency*/, 
	 	    float *ff /*trapezoid key point frequency*/)
/*< initialize with data length >*/
{
   int i1,i2,i3,i4;	
   /* determine frequency sampling (for real to complex FFT) */
   nfft = 2*kiss_fft_next_fast_size((n+1)/2);
   nw  = nfft/2+1;
   dw  = 2.*SF_PI/nfft;
   cdata = (kiss_fft_cpx*) sf_complexalloc(nw);
   shape = sf_floatalloc(nw);
   forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
   invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
   if(NULL==forw || NULL==invs)
	sf_error("%s: KISS FFT allocation problem", __FILE__);

   i1=(ff[0]/fn)*(nw-1)+1; 
   i2=(ff[1]/fn)*(nw-1)+1;
   i3=(ff[2]/fn)*(nw-1)+1; 
   i4=(ff[3]/fn)*(nw-1)+1;
 
   for(iw=0;iw<nw;iw++)
   {
	if(iw<i1) shape[iw]=0;
	if(iw>=i1 && iw<i2) shape[iw]=(iw-i1)/(i2-i1)/nfft;
	if(iw>=i2 && iw<i3) shape[iw]=1/nfft;
	if(iw>=i3 && iw<i4) shape[iw]=(i4-iw)/(i4-i3)/nfft;
	if(iw>=i4) shape[iw]=0;
   }
}

void trapezoid_apply(float *xx /*input time series*/,
		     float *yy /*output time seires*/)
/*< applying trapezoid bandpass filter >*/
{
    kiss_fftr(forw, xx, cdata);
    for(iw=0;iw<nw;iw++)
	cdata[iw] = sf_crmul(cdata[iw],shape[iw]);
    kiss_fftri(invs,cdata,yy);

}

