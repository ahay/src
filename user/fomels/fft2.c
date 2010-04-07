/* 2-D FFT interface */
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

static int n1, n2, nk;
static float wt;
static kiss_fftr_cfg cfg, icfg;
static kiss_fft_cfg cfg2, icfg2;
static kiss_fft_cpx **tmp, *ctrace;
static sf_complex *trace;

int fft2_init(int nx,   int ny   /* input data size */, 
	      int *nx2, int *ny2 /* padded data size */)
/*< initialize >*/
{
    int i2;

    nk = kiss_fft_next_fast_size((nx+1)/2)+1;
    n1 = 2*(nk-1);

    cfg  = kiss_fftr_alloc(n1,0,NULL,NULL);
    icfg = kiss_fftr_alloc(n1,1,NULL,NULL);

    *nx2 = n1;

    n2 = kiss_fft_next_fast_size(ny);

    cfg2  = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(n2,1,NULL,NULL);

    *ny2 = n2;

    tmp = (kiss_fft_cpx **) sf_alloc(n2,sizeof(*tmp));
    tmp[0] = (kiss_fft_cpx *) sf_alloc(nk*n2,sizeof(kiss_fft_cpx));
    for (i2=0; i2 < n2; i2++) {
	tmp[i2] = tmp[0]+i2*nk;
    }

    trace = sf_complexalloc(n2);
    ctrace = (kiss_fft_cpx *) trace;

    wt =  1.0/(n2*n1);

    return (nk*n2);
}

void fft2(float *inp      /* [n1*n2] */, 
	  sf_complex *out /* [nk*n2] */)
/*< 2-D FFT >*/
{
    int i1, i2;

    for (i2=0; i2 < n2; i2++) {
	kiss_fftr (cfg,inp+i2*n1,tmp[i2]);
    }

    /* FFT centering */
    for (i2=1; i2<n2; i2+=2) {
	for (i1=0; i1<nk; i1++) {
	    tmp[i2][i1] = sf_cneg(tmp[i2][i1]);
	}
    }

    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(cfg2,tmp[0]+i1,ctrace,nk);
		
	/* Transpose */
	for (i2=0; i2<n2; i2++) {
	    out[i2*nk+i1] = trace[i2];
	}
    }
}


void ifft2(float *out     /* [n1*n2] */, 
	  sf_complex *inp /* [nk*n2] */)
/*< 2-D inverse FFT >*/
{
    int i1, i2;

    for (i1=0; i1 < nk; i1++) {
	kiss_fft_stride(icfg2,(kiss_fft_cpx *) (inp+i1),ctrace,nk);
		
	for (i2=0; i2<n2; i2++) {
	    tmp[i2][i1] = sf_crmul(ctrace[i2],i2%2? -wt: wt);
	}
    }
    for (i2=0; i2 < n2; i2++) {
	kiss_fftri(icfg,tmp[i2],out+i2*n1);
    }
}
