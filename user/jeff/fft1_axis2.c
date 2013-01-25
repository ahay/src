/* 1-D FFT interface with FFT/IFFT over one axis */
/*
  Copyright (C) 2012 The University of Western Australia
  
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
#include "fft1_axis2.h"
/*^*/

static int nx,nk;
//static sf_complex *trace1, *trace2;
static kiss_fft_cfg cfg1, icfg1;
//static kiss_fft_cpx *ctrace1, *ctrace2;

static float fftscale;

void fft1_axis2_init(int nx_in,int nk_in) 
/*< initialize >*/ 
{
	nx = nx_in;
	nk = nk_in;
    cfg1  = kiss_fft_alloc(nx,0,NULL,NULL);
    icfg1 = kiss_fft_alloc(nk,1,NULL,NULL);
    fftscale = 1./nx;
    
//    trace1 = sf_complexalloc(nx);
//    trace2 = sf_complexalloc(nk);
//    ctrace1 = (kiss_fft_cpx *) trace1;
//    ctrace2 = (kiss_fft_cpx *) trace2;
}

//void fft1_axis2(sf_complex **dat, int id)
void fft1_axis2(kiss_fft_cpx **dat, int id)
/*< 1D FORWARD FOURIER TRANSFORM OVER ONE AXIS >*/
{
    int ix;

	/* PULL OUT AT PROPER PLACE */
//	for (int ix=0; ix < nx; ix++) {
//		trace1[ix] = dat[id][ix];
//	}

	/* FFT on single trace */
    kiss_fft(cfg1,dat[id],dat[id]);
    
    
//    kiss_fft_stride (cfg1,ctrace1,ctrace2,1);

	for (ix=0; ix < nx; ix++) {
	    dat[id][ix] = sf_crmul(dat[id][ix],fftscale);
	}

	/* RETURN BACK TO PROPER PLACE */
//    for (int ix=0; ix < nk; ix++) {
//    	dat[id][ix] = trace2[ix];
//    }
}

//void ifft1_axis2(sf_complex **dat ,int id)
void ifft1_axis2(kiss_fft_cpx **dat ,int id)
/*< 1D INVERSE FOURIER TRANSFORM OVER ONE AXIS >*/
{
	/* PULL OUT AT PROPER PLACE */
//	for (int ix=0; ix < nx; ix++) {
//		trace2[ix] = dat[id][ix];
//	}

	/* IFFT on single trace */	
	kiss_fft (icfg1,dat[id],dat[id]);
//	kiss_fft_stride (icfg1,ctrace2,ctrace1,1);

	/* RETURN BACK TO PROPER PLACE */
//    for (int ix=0; ix < nx; ix++) {
//    	dat[id][ix] = trace1[ix];
//    }
}
