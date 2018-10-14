/* Chain of symmetric Fourier weighting and scaling */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

static int n, nf, nw;
static float *w, *wf, *tmp1, *tmp2, *x1, *x2, *s;
static kiss_fftr_cfg forw, invs;
static kiss_fft_cpx *cdata1, *cdata2;

void sfchain_init(int    n1  /* data size */, 
		  int    nw1 /* Fourier data size */,
		  float *w1  /* weight [n] */,
		  float *wf1 /* Fourier weight [nw] */,
		  float *y1  /* intermediate [n] */,
		  float *y2  /* intermediate [n] */,
		  float *src /* source [n] */)
/*< initialize >*/
{
    n = n1;
    nw = nw1;
    w = w1;
    wf = wf1;
    x1 = y1;
    x2 = y2;
    s = src;

    nf = 2*(nw-1);

    forw = kiss_fftr_alloc(nf,0,NULL,NULL);
    invs = kiss_fftr_alloc(nf,1,NULL,NULL);

    tmp1 = sf_floatalloc(nf);    
    tmp2 = sf_floatalloc(nf);  

    cdata1 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
    cdata2 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
}

void sfchain_close(void)
/*< clean allocated storage >*/
{
    free(tmp1);
    free(tmp2);
    free(cdata1);
    free(cdata2);
    free(forw);
    free(invs);
}

void sfchain_res(const float *t /* target */,
		 float *r       /* residual */)
/*< apply the chain operator >*/
{
    int i, iw;

    for (i=0; i < n; i++) {
	tmp2[i] = x1[i];
    }
    for (i=n; i < nf; i++) {
	tmp2[i] = 0.0f;
    }
    kiss_fftr(forw, tmp2, cdata2);
    for (iw=0; iw < nw; iw++) {
	cdata2[iw]=sf_crmul(cdata2[iw],wf[iw]/nf);
    }
    kiss_fftri(invs, cdata2, tmp2);
    for (i=0; i < n; i++) {
	r[i]     = x1[i]-w[i]*s[i];
	r[n+i]   = x2[i]-tmp2[i];
	r[2*n+i] = t[i]-w[i]*x2[i];
    }
}

void sfchain_apply(float *y)
/*< apply the chain operator >*/
{
    int i, iw;

    for (i=0; i < n; i++) {
	tmp2[i] = w[i]*s[i];
    }
    for (i=n; i < nf; i++) {
	tmp2[i] = 0.0f;
    }
    kiss_fftr(forw, tmp2, cdata2);
    for (iw=0; iw < nw; iw++) {
	cdata2[iw]=sf_crmul(cdata2[iw],wf[iw]/nf);
    }
    kiss_fftri(invs, cdata2, tmp2);
    for (i=0; i < n; i++) {
	y[i] = w[i]*tmp2[i];
    }
}


void sfchain_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear operator >*/
{
    int i, iw;

    if (nx != 3*n+nw || ny != 3*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);
    
    if (adj) {
	for (i=0; i < n; i++) {
	    tmp1[i] = y[n+i];
	    tmp2[i] = x1[i];
	    x[n+i] += w[i]*y[2*n+i] - y[n+i];
	    x[2*n+i] += x2[i]*y[2*n+i];
	}
	for (i=n; i < nf; i++) {
	    tmp1[i] = 0.0f;
	    tmp2[i] = 0.0f;
	}
	kiss_fftr(forw, tmp1, cdata1);
	kiss_fftr(forw, tmp2, cdata2);
	for (iw=0; iw < nw; iw++) {
	    x[3*n+iw] += sf_crealf(sf_cmul(cdata1[iw],sf_conjf(cdata2[iw])))/nf;
	    cdata1[iw]=sf_crmul(cdata1[iw],wf[iw]/nf);
	}
	kiss_fftri(invs, cdata1, tmp1);
	for (i=0; i < n; i++) {
	    x[2*n+i] += s[i]*y[i];
	    x[i] += tmp1[i] - y[i];
	}
    } else {
	for (i=0; i < n; i++) {
	    y[i] += s[i]*x[2*n+i] - x[i];
	    tmp1[i] = x[i];
	    tmp2[i] = x1[i];
	}
	for (i=n; i < nf; i++) {
	    tmp1[i] = 0.0f;
	    tmp2[i] = 0.0f;
	}
	kiss_fftr(forw, tmp1, cdata1);
	kiss_fftr(forw, tmp2, cdata2);
	for (iw=0; iw < nw; iw++) {
	    cdata1[iw]=sf_crmul(cdata1[iw],wf[iw]/nf);
	    cdata2[iw]=sf_crmul(cdata2[iw],x[3*n+iw]/nf);
	}
	kiss_fftri(invs, cdata1, tmp1);
	kiss_fftri(invs, cdata2, tmp2);
	for (i=0; i < n; i++) {
	    y[n+i] += tmp1[i] + tmp2[i] - x[n+i];
	    y[2*n+i] += x2[i]*x[2*n+i] + w[i]*x[n+i];
	}
    }
}
	
	
