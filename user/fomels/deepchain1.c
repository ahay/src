/* Deep symmetric chain with 1D Fourier*/
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


static int nf, nw, nt, nx, ntx, nwx;
static float *w, *wf, *x1, *x2, *x3, *x4, *s;
static float **tmp1, **tmp2, **tmp3, **tmp4;

static kiss_fftr_cfg forw, invs;
static kiss_fft_cpx *cdata1, *cdata2, *cdata3, *cdata4;

void sfdeepchain1_init(int    n1  /* trace length */,
		  int    n2  /* number of traces */,
		  int    nw1 /* Fourier trace length */,
		  float *w1  /* space weight [n1*n2] */,
		  float *wf1 /* Fourier weight [nw*n2] */,
		  float *y1  /* intermediate [n1*n2] */,
		  float *y2  /* intermediate [n1*n2] */,
		  float *y3  /* intermediate [n1*n2] */,
		  float *y4  /* intermediate [n1*n2] */,	
		  float *src /* source [n1*n2] */)
/*< initialize >*/
{
	nt = n1;
	nx = n2;
	nw = nw1;
	w = w1;
	wf = wf1;
	x1 = y1;
	x2 = y2;
	x3 = y3;
	x4 = y4;
	s = src;

	nf = 2*(nw1-1);
	ntx = nt*nx;
	nwx = nw*nx;

    forw = kiss_fftr_alloc(nf,0,NULL,NULL);
    invs = kiss_fftr_alloc(nf,1,NULL,NULL);

    tmp1 = sf_floatalloc2(nf,nx);    
    tmp2 = sf_floatalloc2(nf,nx);  
    tmp3 = sf_floatalloc2(nf,nx);  
    tmp4 = sf_floatalloc2(nf,nx);  

    cdata1 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
    cdata2 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
    cdata3 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
    cdata4 =  (kiss_fft_cpx*) sf_alloc(nw,sizeof(kiss_fft_cpx));
}

void sfdeepchain1_close(void)
/*< clean allocated storage >*/
{
    free(*tmp1);
    free(*tmp2);
    free(*tmp3);
    free(*tmp4);
    free(tmp1);
    free(tmp2);
    free(tmp3);
    free(tmp4);
    free(cdata1);
    free(cdata2);
    free(cdata3);
    free(cdata4);
    free(forw);
    free(invs);
}

void sfdeepchain1_res(const float *t /* target */,
		 float *r       /* residual */)
/*< apply the chain operator >*/
{
    int i, iw, i2, i1;

    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    tmp1[i2][i1] = x1[i1+i2*nt];
	    tmp2[i2][i1] = x3[i1+i2*nt];
	}
	for (i1=nt; i1 < nf; i1++) {
	    tmp1[i2][i1] = 0.0f;
	    tmp2[i2][i1] = 0.0f;
	}
	kiss_fftr(forw, tmp1[i2], cdata1);
	kiss_fftr(forw, tmp2[i2], cdata2);
	for (iw=0; iw < nw; iw++) {
	    cdata1[iw]=sf_crmul(cdata1[iw],wf[iw+i2*nw]/nf);
	    cdata2[iw]=sf_crmul(cdata2[iw],wf[iw+i2*nw]/nf);
	}
	kiss_fftri(invs, cdata1, tmp1[i2]);
	kiss_fftri(invs, cdata2, tmp2[i2]);
	} /* End of i2 iteration */

    for (i=0; i < ntx; i++) {
	i1 = i%nt;
	i2 = i/nt;
	
	r[i]     = x1[i]-w[i]*s[i];
	r[ntx+i]   = x2[i]- tmp1[i2][i1];
	r[2*ntx+i] = x3[i] - w[i]*x2[i];
	r[3*ntx+i] = x4[i] - tmp2[i2][i1];
	r[4*ntx+i] = t[i] - w[i]*x4[i]; 
    }
}


void sfdeepchain1_apply(float *y)
/*< apply the chain operator >*/
{
    int i2, i1, i, iw;

    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    i = i1+i2*nt;
	    tmp1[i2][i1] = w[i]*s[i];
	}
	for (i1=ntx; i1 < nf; i1++) {
	    tmp1[i2][i1] = 0.0f;
	}
	kiss_fftr(forw, tmp1[i2], cdata1);

	for (iw=0; iw < nw; iw++) {
	    cdata1[iw]=sf_crmul(cdata1[iw],wf[iw+i2*nw]/nf);
	}

	kiss_fftri(invs, cdata1, tmp1[i2]);

	for (i1=0; i1 < nt; i1++) {
	    i = i1+i2*nt;
	    tmp1[i2][i1] = w[i]*tmp1[i2][i1];
	}
	kiss_fftr(forw, tmp1[i2], cdata1);

	for (iw=0; iw < nw; iw++) {
	    cdata1[iw]=sf_crmul(cdata1[iw],wf[iw+i2*nw]/nf);
	}

	kiss_fftri(invs, cdata1, tmp1[i2]);
	} /* End of i2 iteration */

    for (i=0; i < ntx; i++) {
	i1 = i%nt;
	i2 = i/nt;
	
	y[i] = w[i]*tmp1[i2][i1];
    }
}


void sfdeepchain1_deconimg(const float *t ,float *lsmig, float* spaceW, float* freqW)
/*<apply the chain to the target (first mig) equavalent to ls migration >*/
{
    int i2, i1, i, iw;

    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    i = i1+i2*nt;
	    tmp1[i2][i1] = spaceW[i]*t[i];
	}
	for (i1=ntx; i1 < nf; i1++) {
	    tmp1[i2][i1] = 0.0f;
	}
	kiss_fftr(forw, tmp1[i2], cdata1);

	for (iw=0; iw < nw; iw++) {
	    cdata1[iw]=sf_crmul(cdata1[iw],freqW[iw+i2*nw]/nf);
	}

	kiss_fftri(invs, cdata1, tmp1[i2]);


	for (i1=0; i1 < nt; i1++) {
	    i = i1+i2*nt;
	    tmp1[i2][i1] = spaceW[i]*tmp1[i2][i1];
	}

	kiss_fftr(forw, tmp1[i2], cdata1);

	for (iw=0; iw < nw; iw++) {
	    cdata1[iw]=sf_crmul(cdata1[iw],freqW[iw+i2*nw]/nf);
	}

	kiss_fftri(invs, cdata1, tmp1[i2]);
	} /* End of i2 iteration */

    for (i=0; i < ntx; i++) {
	i1 = i%nt;
	i2 = i/nt;
	
	lsmig[i] = spaceW[i]*tmp1[i2][i1];
    }
}

void sfdeepchain1_lop (bool adj, bool add, int nxx, int nyy, float* x, float* y) 
/*< linear operator >*/
{

    int i, iw, i1, i2;
    float buf1, buf2, scale;


    if (nxx != 5*ntx+nwx || nyy != 5*ntx) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nxx,nyy,x,y);
    
    if (adj) { /* Adjoint */

    /* Initialization */ 
	for (i=0; i < ntx; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    tmp1[i2][i1] = y[ntx+i]; 
	    tmp2[i2][i1] = y[3*ntx+i]; 
	    tmp3[i2][i1] = x1[i];
	    tmp4[i2][i1] = x3[i];

	}

	for (i2=0; i2 < nx; i2++) {
	    for (i1=nt; i1 < nf; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
		tmp3[i2][i1] = 0.0f;
		tmp4[i2][i1] = 0.0f;
	    }
	    kiss_fftr(forw, tmp1[i2], cdata1);
	    kiss_fftr(forw, tmp2[i2], cdata2);
	    kiss_fftr(forw, tmp3[i2], cdata3);
	    kiss_fftr(forw, tmp4[i2], cdata4);
	    for (iw=0; iw < nw; iw++) {

        if(iw==0 || iw == nw-1){scale = 1.0;}
        else{scale = 2.0;}

	    buf1 = sf_crealf(sf_cmul(cdata1[iw],sf_conjf(cdata3[iw])))/nf;
	    buf2 = sf_crealf(sf_cmul(cdata2[iw],sf_conjf(cdata4[iw])))/nf;

		x[5*ntx+iw+i2*nw] += scale*(buf1 + buf2);


		cdata1[iw]=sf_crmul(cdata1[iw],wf[iw+i2*nw]/nf);
		cdata2[iw]=sf_crmul(cdata2[iw],wf[iw+i2*nw]/nf);
		}
	    kiss_fftri(invs, cdata1, tmp1[i2]);
	    kiss_fftri(invs, cdata2, tmp2[i2]);


	} /* End of i2 iterations */

	for (i=0; i < ntx; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    x[i] += tmp1[i2][i1] - y[i];
	    x[ntx+i] += w[i]*y[2*ntx+i] - y[ntx+i];
	    x[2*ntx+i] += tmp2[i2][i1] - y[2*ntx+i];
	    x[3*ntx+i] += w[i]*y[4*ntx+i] - y[3*ntx+i];
	    x[4*ntx+i] += s[i]*y[i] + x2[i]*y[2*ntx+i] + x4[i]*y[4*ntx+i];

	}


    } else { /* Forward */

    /* Initialization */ 
	for (i=0; i < ntx; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    tmp1[i2][i1] = x[i];
	    tmp2[i2][i1] = x1[i];
	    tmp3[i2][i1] = x[2*ntx+i];
	    tmp4[i2][i1] = x3[i];

	}

	for (i2=0; i2 < nx; i2++) {

	    for (i1=nt; i1 < nf; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
		tmp3[i2][i1] = 0.0f;
		tmp4[i2][i1] = 0.0f;
	    }

	    kiss_fftr(forw, tmp1[i2], cdata1);
	    kiss_fftr(forw, tmp2[i2], cdata2);
	    kiss_fftr(forw, tmp3[i2], cdata3);
	    kiss_fftr(forw, tmp4[i2], cdata4);

	    for (iw=0; iw < nw; iw++) {

		cdata1[iw]=sf_crmul(cdata1[iw],wf[iw+i2*nw]/nf);
		cdata2[iw]=sf_crmul(cdata2[iw],x[5*ntx+iw+i2*nw]/nf);
		cdata3[iw]=sf_crmul(cdata3[iw],wf[iw+i2*nw]/nf);
		cdata4[iw]=sf_crmul(cdata4[iw],x[5*ntx+iw+i2*nw]/nf);
		}
	    kiss_fftri(invs, cdata1, tmp1[i2]);
	    kiss_fftri(invs, cdata2, tmp2[i2]);
	    kiss_fftri(invs, cdata3, tmp3[i2]);
	    kiss_fftri(invs, cdata4, tmp4[i2]);


	} /* End of i2 iterations */

	for (i=0; i < ntx; i++) {
	    i1 = i%nt;
	    i2 = i/nt;

	    y[i] += s[i]*x[4*ntx+i] - x[i];
	    y[ntx+i] += tmp1[i2][i1] + tmp2[i2][i1] - x[ntx+i] ;

	    y[2*ntx+i] += w[i]*x[ntx+i] + x2[i]*x[4*ntx+i] - x[2*ntx+i];
	    y[3*ntx+i] += tmp3[i2][i1] + tmp4[i2][i1] - x[3*ntx+i];

	    y[4*ntx+i] += w[i]*x[3*ntx+i] + x4[i]*x[4*ntx+i];
	}
	} /* End of else */
}

