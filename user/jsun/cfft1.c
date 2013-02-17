/* Complex 1-D FFT interface, complex to complex fft and ifft */
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

static int n1, nk;
static float wt;

static sf_complex *cc=NULL;

static kiss_fft_cfg cfg1, icfg1;
static kiss_fft_cpx *tmp;
static sf_complex *tmp2;

int fft1_init(int pad1           /* padding on the first axis */,
	      int nx             /* input data size */, 
	      int *nx2           /* padded data size */)
/*< initialize >*/
{
	

    nk = n1 = kiss_fft_next_fast_size(nx*pad1);
    
    cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
    icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);
    
    cc = sf_complexalloc(n1);
     	
    tmp2 = sf_complexalloc(nk);
    tmp  = (kiss_fft_cpx *) tmp2;

    *nx2 = n1;
	
    wt =  1.0/n1;
	
    return (nk);
}

void fft1(sf_complex *inp /* [n1] */, 
	  sf_complex *out /* [nk] */)
/*< 1-D FFT >*/
{
    int i1;
//    kiss_fft_cpx *inp1;
//    inp1 = (kiss_fft_cpx *) inp;

    /* FFT centering */
    for (i1=0; i1<n1; i1++) {

#ifdef SF_HAS_COMPLEX_H
	cc[i1] = i1%2? inp[i1]:(-1*inp[i1]);
#else
	cc[i1] = i1%2? inp[i1]:sf_cneg(inp[i1]);
#endif
    }
    

    kiss_fft_stride(cfg1,(kiss_fft_cpx *) cc,tmp,1);
	
    for (i1=0; i1<n1; i1++) {
	    out[i1] = tmp2[i1];
	}
}

void ifft1_allocate(sf_complex *inp /* [nk*n2] */)
/*< allocate inverse transform >*/
{
#ifdef SF_HAS_FFTW

#endif
}

void ifft1(sf_complex *out /* [n1] */, 
	   sf_complex *inp /* [nk] */)
/*< 1-D inverse FFT >*/
{
    int i1;

    kiss_fft_stride(icfg1,(kiss_fft_cpx *) inp,(kiss_fft_cpx *) cc,1);
    
    /* FFT centering and normalization */
    for (i1=0; i1<n1; i1++) {

#ifdef SF_HAS_COMPLEX_H
	out[i1] = (i1%2? wt:-wt) * cc[i1];
#else
	out[i1] = sf_crmul(cc[i1],(i1%2? wt:-wt));
#endif
    }
}
