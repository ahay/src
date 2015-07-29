/* S transform */
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
#include <stdio.h>
#include <math.h>
#include "st.h"

float gauss(int n, int m);

void st (int len            /* data size              */, 
	 float d1           /* data sampling          */, 
	 int lo             /* low frequency          */, 
	 int hi             /* high frequency         */, 
	 float *data        /* input [len]            */, 
	 sf_complex *result /* output [len*(hi-lo+1)] */)
/*< Forward S transform >*/
{
    int i, i1, k, l2, nw;
    float s, *g;

    kiss_fft_cpx *d, *pp, *qq;
    kiss_fft_cfg tfft, itfft;

    nw = 2*kiss_fft_next_fast_size((len+1)/2);
    tfft = kiss_fft_alloc(nw,0,NULL,NULL);
    itfft = kiss_fft_alloc(nw,1,NULL,NULL);

    pp = (kiss_fft_cpx*) sf_complexalloc(nw);
    qq = (kiss_fft_cpx*) sf_complexalloc(nw);
    d =  (kiss_fft_cpx*) sf_complexalloc(nw);
    g = sf_floatalloc(nw);

    s = 0.;
    for (i = 0; i < len; i++) {
		d[i].r = data[i];
		d[i].i = 0.;
		s += data[i];
	}
    s /= len;
    
    for (i=len; i < nw; i++) {
		d[i].r = 0.;
		d[i].i = 0.;
    }		    

    kiss_fft_stride (tfft,d,pp,1);
    
    l2 = (nw+1)/2;
    for (i=1; i < l2; i++) {
		pp[i].r *= 2.;
		pp[i].i *= 2.;
    }
    l2 = nw/2+1;
    for (i=l2; i < nw; i++) {
		pp[i].r = 0.;
		pp[i].i = 0.;
    }

    for (i1=lo; i1 <= hi; i1++) {
		if (0 == i1) {
			for (i=0; i < len; i++) {
				result[(i1-lo)*len+i] = sf_cmplx(s,0.);
			}	
		} else {
			g[0] = gauss(i1, 0);
			l2 = nw/2 + 1;
			for (i=1; i < l2; i++) {
				g[i] = g[nw-i] = gauss(i1, i);
			}
	    
			for (i=0; i < nw; i++) {
				s = g[i];
				k = i1 + i;
				if (k >= nw) k -= nw;
				qq[i].r = pp[k].r * s;
				qq[i].i = pp[k].i * s;
			}
	    
			kiss_fft_stride(itfft,qq,d,1);
	    
			for (i=0; i < len; i++) {
				result[(i1-lo)*len+i] = sf_cmplx(d[i].r/len,d[i].i/len);
			}
		}	    
    }
    free(pp);
    free(qq);
    free(d);
    free(g);
}

void ist (int len            /* data size              */, 
	  float d1           /* data sampling          */, 
	  int lo             /* low frequency          */, 
	  int hi             /* high frequency         */, 
	  float *result      /* output [len]           */, 
	  sf_complex *data   /* input [len*(hi-lo+1)]  */)
/*< Inverse S transform >*/
{
    int i, i1, l2, nw;

    kiss_fft_cpx *d, *pp;
    kiss_fft_cfg itfft;

    nw = 2*kiss_fft_next_fast_size((len+1)/2);    
    itfft = kiss_fft_alloc(nw,1,NULL,NULL);

    pp = (kiss_fft_cpx*) sf_complexalloc(nw);
    d =  (kiss_fft_cpx*) sf_complexalloc(nw);

    for (i=0; i < nw; i++) {
		pp[i].r = 0.;
		pp[i].i = 0.;
    }

    for (i1=lo; i1 <= hi; i1++) {
		for (i=0; i < len; i++) {
			pp[i1-lo].r += crealf(data[(i1-lo)*len+i]);
			pp[i1-lo].i += cimagf(data[(i1-lo)*len+i]);
		}
    }
 
    l2 = (nw+1)/2;
    for (i=1; i < l2; i++) {
		pp[i].r /= 2.;
		pp[i].i /= 2.;
    }
    l2 = nw/2+1;
    for (i=l2; i < nw; i++) {
		pp[i].r = pp[nw-i].r;
		pp[i].i = -pp[nw-i].i;
    }
    kiss_fft_stride(itfft,pp,d,1);
	    
    for (i=0; i < len; i++) {
		result[i] = d[i].r/len;
    }
    free(pp);
    free(d);
}

float gauss(int n, int m)
/*< Fourier Transform of a Gaussian >*/
{
    return exp(-2.*SF_PI*SF_PI*m*m/(n*n));
}



/* 	$Id: st.c 8858 2012-07-23 16:33:06Z saragiotis $	 */

