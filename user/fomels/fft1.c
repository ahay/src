/* 1-D FFT interface */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

static kiss_fftr_cfg cfg, icfg;

int fft1_init(int n1  /* input data size */, 
	      int *n2 /* padded data size */)
/*< initialize >*/
{
    int nk;

    nk = kiss_fft_next_fast_size((n1+1)/2)+1;
    n1 = 2*(nk-1);

    cfg  = kiss_fftr_alloc(n1,0,NULL,NULL);
    icfg = kiss_fftr_alloc(n1,1,NULL,NULL);

    *n2 = n1;
    return nk;
}

void fft1(float *inp      /* [n2] */, 
	 sf_complex *out /* [nk] */)
/*< 1-D FFT >*/
{
    kiss_fftr (cfg,inp,(kiss_fft_cpx *) out);
}


void ifft1(float *out      /* [n2] */, 
	  sf_complex *inp /* [nk] */)
/*< 1-D inverse FFT >*/
{
    kiss_fftri(icfg,(kiss_fft_cpx *) inp,out);
}
