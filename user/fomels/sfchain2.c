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

static int n, nf, nw, nt, nx;
static float *w, *wf, **tmp1, **tmp2, *x1, *x2, *s;
static kiss_fftr_cfg forw, invs;
static kiss_fft_cpx *cdata1, *cdata2;

void sfchain2_init(int    n1  /* trace length */,
		  int    n2  /* number of traces */,
		  int    nw1 /* Fourier trace length */,
		  float *w1  /* weight [n1*n2] */,
		  float *wf1 /* Fourier weight [nw*n2] */,
		  float *y1  /* intermediate [n1*n2] */,
		  float *y2  /* intermediate [n1*n2] */,
		  float *src /* source [n1*n2] */)
/*< initialize >*/
{
    nt = n1;
    nx = n2;
    
    n = n1*n2;
    nw = nw1;
    w = w1;
    wf = wf1;
    x1 = y1;
    x2 = y2;
    s = src;

    nf = 2*(nw1-1);

    forw = kiss_fftr_alloc(nf,0,NULL,NULL);
    invs = kiss_fftr_alloc(nf,1,NULL,NULL);

    tmp1 = sf_floatalloc2(nf,nx);    
    tmp2 = sf_floatalloc2(nf,nx);  

    cdata1 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
    cdata2 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
}

void sfchain2_close(void)
/*< clean allocated storage >*/
{
    free(*tmp1);
    free(*tmp2);
    free(tmp1);
    free(tmp2);
    free(cdata1);
    free(cdata2);
    free(forw);
    free(invs);
}

void sfchain2_res(const float *t /* target */,
		 float *r       /* residual */)
/*< apply the chain operator >*/
{
    int i, iw, i2, i1;

    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    tmp2[i2][i1] = x1[i1+i2*nt];
	}
	for (i1=nt; i1 < nf; i1++) {
	    tmp2[i2][i1] = 0.0f;
	}
	kiss_fftr(forw, tmp2[i2], cdata2);
	for (iw=0; iw < nw; iw++) {
	    cdata2[iw]=sf_crmul(cdata2[iw],wf[iw+i2*nw]/nf);
	}
	kiss_fftri(invs, cdata2, tmp2[i2]);
    }
    for (i=0; i < n; i++) {
	i1 = i%nt;
	i2 = i/nt;
	
	r[i]     = x1[i]-w[i]*s[i];
	r[n+i]   = x2[i]-tmp2[i2][i1];
	r[2*n+i] = t[i]-w[i]*x2[i];
    }
}

void sfchain2_apply(float *y)
/*< apply the chain operator >*/
{
    int i2, i1, i, iw;

    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    i = i1+i2*nt;
	    tmp2[i2][i1] = w[i]*s[i];
	}
	for (i1=n; i1 < nf; i1++) {
	    tmp2[i2][i1] = 0.0f;
	}
	kiss_fftr(forw, tmp2[i2], cdata2);
	for (iw=0; iw < nw; iw++) {
	    cdata2[iw]=sf_crmul(cdata2[iw],wf[iw+i2*nw]/nf);
	}
	kiss_fftri(invs, cdata2, tmp2[i2]);
    }
    for (i=0; i < n; i++) {
	i1 = i%nt;
	i2 = i/nt;
	
	y[i] = w[i]*tmp2[i2][i1];
    }
}


void sfchain2_lop (bool adj, bool add, int nxx, int nyy, float* x, float* y) 
/*< linear operator >*/
{
    int i, iw, i1, i2;

    if (nxx != 3*n+nw*nx || nyy != 3*n) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nxx,nyy,x,y);
    
    if (adj) {
	for (i=0; i < n; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    tmp1[i2][i1] = y[n+i];
	    tmp2[i2][i1] = x1[i];
	    x[n+i] += w[i]*y[2*n+i] - y[n+i];
	    x[2*n+i] += x2[i]*y[2*n+i];
	}
	for (i2=0; i2 < nx; i2++) {
	    for (i1=nt; i1 < nf; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
	    }
	    kiss_fftr(forw, tmp1[i2], cdata1);
	    kiss_fftr(forw, tmp2[i2], cdata2);
	    for (iw=0; iw < nw; iw++) {
		x[3*n+iw+i2*nw] +=
		    sf_crealf(sf_cmul(cdata1[iw],sf_conjf(cdata2[iw])))/nf;
		cdata1[iw]=sf_crmul(cdata1[iw],wf[iw+i2*nw]/nf);
	    }
	    kiss_fftri(invs, cdata1, tmp1[i2]);
	}
	for (i=0; i < n; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    x[2*n+i] += s[i]*y[i];
	    x[i] += tmp1[i2][i1] - y[i];
	}
    } else {
	for (i=0; i < n; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    y[i] += s[i]*x[2*n+i] - x[i];
	    tmp1[i2][i1] = x[i];
	    tmp2[i2][i1] = x1[i];
	}
	for (i2=0; i2 < nx; i2++) {
	    for (i1=nt; i1 < nf; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
	    }
	    kiss_fftr(forw, tmp1[i2], cdata1);
	    kiss_fftr(forw, tmp2[i2], cdata2);
	    for (iw=0; iw < nw; iw++) {
		cdata1[iw]=sf_crmul(cdata1[iw],wf[iw+i2*nw]/nf);
		cdata2[iw]=sf_crmul(cdata2[iw],x[3*n+iw+i2*nw]/nf);
	    }
	    kiss_fftri(invs, cdata1, tmp1[i2]);
	    kiss_fftri(invs, cdata2, tmp2[i2]);
	}
	for (i=0; i < n; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    y[n+i] += tmp1[i2][i1] + tmp2[i2][i1] - x[n+i];
	    y[2*n+i] += x2[i]*x[2*n+i] + w[i]*x[n+i];
	}
    }
}
	
	
