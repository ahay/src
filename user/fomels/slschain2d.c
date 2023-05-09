/* Separable LS Chain of symmetric 2D-Fourier weighting and scaling */
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
#include "slschain2d.h"
#include "fft2.h"
#include "cdivn.h"



static int n, nk, nw, nt, nx, nt1,nt2,nx2;
static float *w, *wf, **tmp1, **tmp2, *x1, *x2, *s;
static sf_complex *ctmp1, *ctmp2; /* for 2D-fft */

void sfslschain2d_init(int  n1     /* trace length */,
		       int  n2     /* number of traces */,
		       int  n1pad  /* dim1 for 2-D fft input */,
		       int  n_fftx /* dim2 for 2-D fft input */,
		       int  n_out_fft     /* dim for 2-D FFT output */,
		       float *w1   /* weight [n1*n2] */,
		       float *wf1  /* Fourier weight [nk = nk*n2 (from fft2)) ] */,
		       float *y1   /* intermediate [n1*n2] */,
		       float *y2   /* intermediate [n1*n2] */,
		       float *src  /* source [n1*n2] */)
/*< initialize >*/
{
    nt = n1;
    nx = n2;
    nt1 = n1pad;
    nx2 = n_fftx;
    nk = n_out_fft;
    nt2 = nk/n_fftx;
    nw = nt2;

    /*Model parameters*/
    n = n1*n2;
    w = w1;
    wf = wf1; 
    x1 = y1; 
    x2 = y2; 
    s = src; 

    tmp1 = sf_floatalloc2(nt1,nx2);  
    tmp2 = sf_floatalloc2(nt1,nx2);    

    /*for 2D fft*/
    ctmp1 = sf_complexalloc(nk);
    ctmp2 = sf_complexalloc(nk);

}

void sfslschain2d_close(void)
/*< clean allocated storage >*/
{
    free(*tmp1);
    free(*tmp2);
    free(tmp1);
    free(tmp2);
    free(ctmp1);
    free(ctmp2);
}

void sfslschain2d_res(const float *t /* target */,
		      float *r       /* residual */)
/*< apply the chain operator >*/

{
    int i, i1, i2, ik;

    /* pad with zeros */
    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    tmp2[i2][i1] = x1[i1+i2*nt];
	}
	for (i1=nt; i1 < nt1; i1++) {
	    tmp2[i2][i1] = 0.0f;
	}
    }
    for (i2=nx; i2 < nx2; i2++) {
	for (i1=0; i1 < nt1; i1++) {
	    tmp2[i2][i1] = 0.0f;
	}
    }

    /* forward FFT */

    fft2_allocate(ctmp2);
    ifft2_allocate(ctmp2);
    fft2(tmp2[0],ctmp2);

    /* frequency weight */
    for (ik=0; ik < nk; ik++) {
	ctmp2[ik] *= wf[ik];
    }
    /* inverse FFT */
    ifft2(tmp2[0],ctmp2);
    /* Compute residual r */
    for(i=0; i<nt*nx; i++){
	i1 = i%nt;
	i2 = i/nt;

	r[i]     = x1[i]-w[i]*s[i];
	r[n+i]   = x2[i]-tmp2[i2][i1];
	r[2*n+i] = t[i]-w[i]*x2[i];
    }   
}

void sfslschain2d_apply(float *y)
/*< apply the chain operator >*/

{
    int i, ik, i1, i2;

    /* pad with zeros */
    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    i = i1+i2*nt;
	    tmp2[i2][i1] = w[i]*s[i];
	}
	for (i1=nt; i1 < nt1; i1++) {
	    tmp2[i2][i1] = 0.0f;
	}
    }
    for (i2=nx; i2 < nx2; i2++) {
	for (i1=0; i1 < nt1; i1++) {
	    tmp2[i2][i1] = 0.0f;
	}
    }

    /* forward FFT */
    fft2_allocate(ctmp2);
    fft2(tmp2[0],ctmp2);

    /* frequency weight */
    for (ik=0; ik < nk; ik++) {
	ctmp2[ik] *= wf[ik];
    }
    /* inverse FFT */
    ifft2_allocate(ctmp2);
    ifft2(tmp2[0],ctmp2);

    for (i=0; i < n; i++) {
	i1 = i%nt;
	i2 = i/nt;
	
	y[i] = w[i]*tmp2[i2][i1];
    }
}


void sfslschain2d_lop (bool adj, bool add, int nxx, int nyy, float* x, float* y) 
/*< linear operator >*/
{
    int i, ik, i1, i2;
    int *dc_id, *nyq_id;
    int id;
    float scale;
    dc_id = sf_intalloc(nx);
    nyq_id = sf_intalloc(nx);


    if (nxx != 3*n+nk || nyy != 3*n) sf_error("%s: Wrong size",__FILE__);
    
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

	/* pad with zeros */
	for (i2=0; i2 < nx; i2++) {
	    for (i1=nt; i1 < nt1; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
	    }
	}
	for (i2=nx; i2 < nx2; i2++) {
	    for (i1=0; i1 < nt1; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
	    }
	}

	fft2_allocate(ctmp2);
	fft2(tmp2[0],ctmp2);
	fft2_allocate(ctmp1);
	ifft2_allocate(ctmp1);
	fft2(tmp1[0],ctmp1);

	for(id = 0; id < nx; id++){
	    dc_id[id] = id*nw;
	    nyq_id[id] = (id+1)*nw - 1;
	}

	for (ik=0; ik < nk; ik++) {
	    scale = 2.0;
	    for(id=0; id<nx; id++)
	    { 
		if(ik==dc_id[id] || ik == nyq_id[id]){
		    scale = 1.0;
		    break;
		}
	    }

	    x[3*n+ik] += 
		scale*crealf(ctmp1[ik]*conjf(ctmp2[ik])/(nt1*nx2));

	    ctmp1[ik] *= wf[ik];
	}	
	ifft2(tmp1[0],ctmp1);

	for (i=0; i < n; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    x[2*n+i] += s[i]*y[i];
	    x[i] += tmp1[i2][i1] - y[i];
	}
    } else { /*Forward*/

	for (i=0; i < n; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    
	    y[i] += s[i]*x[2*n+i] - x[i];
	    tmp1[i2][i1] = x[i];
	    tmp2[i2][i1] = x1[i];
	}
	/* pad with zeros */
	for (i2=0; i2 < nx; i2++) {
	    for (i1=nt; i1 < nt1; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
	    }
	}
	for (i2=nx; i2 < nx2; i2++) {
	    for (i1=0; i1 < nt1; i1++) {
		tmp1[i2][i1] = 0.0f;
		tmp2[i2][i1] = 0.0f;
	    }
	}

	fft2_allocate(ctmp1);	
	fft2(tmp1[0],ctmp1);	
	for(ik=0; ik<nk;ik++){
	    ctmp1[ik] *=wf[ik];
	}
	ifft2(tmp1[0],ctmp1);
	fft2_allocate(ctmp2);
	ifft2_allocate(ctmp2);
	fft2(tmp2[0],ctmp2);
	for(ik=0; ik<nk;ik++){
	    ctmp2[ik] *=x[3*n+ik];
	}
	ifft2(tmp2[0],ctmp2);

	for (i=0; i < n; i++) {
	    i1 = i%nt;
	    i2 = i/nt;
	    y[n+i] += tmp1[i2][i1] + tmp2[i2][i1] - x[n+i];
	    y[2*n+i] += x2[i]*x[2*n+i] + w[i]*x[n+i];
	}
    }
}

void sfslschain2d_wf(float* w_freq,float* tgt, float *src, float *w_time, int frect1, int frect2)
/*<Solve for Wf using Separable LS approach >*/
{
    int i, i1,i2, id, dim, nd, niter, rect[2], ndim[SF_MAX_DIM];
    float norm, a, **num, **den;
    sf_complex *c_num, *c_den, *cw_freq;


    num = sf_floatalloc2(nt1,nx2);  
    den = sf_floatalloc2(nt1,nx2);    

    /*for 2D fft*/
    c_num = sf_complexalloc(nk);
    c_den = sf_complexalloc(nk);
    cw_freq = sf_complexalloc(nk);

    /* initialize t/w and w*s */ 
    for (i2=0; i2 < nx; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    i = i1+i2*nt;
	    num[i2][i1] = tgt[i]/w_time[i];
	    den[i2][i1] = w_time[i]*src[i];
	}
	for (i1=nt; i1 < nt1; i1++) {
	    num[i2][i1] = 0.0f;
	    den[i2][i1] = 0.0f;
	}
    }
    for (i2=nx; i2 < nx2; i2++) {
	for (i1=0; i1 < nt1; i1++) {
	    num[i2][i1] = 0.0f;
	    den[i2][i1] = 0.0f;
	}
    }

    /* FFT of num and den */
    fft2_allocate(c_num);	
    fft2(num[0],c_num);

    fft2_allocate(c_den);
    fft2(den[0],c_den);

//////////////////////////////////////////
    /* Complex smooth division*/
    /* Imitate Mcdivn.c */
    niter = 100;
    //dim = sf_filedims (fnum,n);
    dim = 2; // output of 2d fft?
    nd = nk; //data size?
    rect[0] = frect1; //previously used for chain lop
    rect[1] = frect2;
    // Not sure still about ndim (n in Mcdivn.c)
    ndim[0] = nk;
    ndim[1] = 1;

/*
// how to initialize smoothing radius here ?? //
for (i=0; i < dim; i++) {
snprintf(key,6,"rect%d",i+1);
if (!sf_getint(key,rect+i)) rect[i]=1;
//( rect#=(1,1,...) smoothing radius on #-th axis )
//nd *= n[i]; n1*n2*... but here we have nk(=nk*n2) already
} */
////////////////////////////////////////

    cdivn_init(dim, nd, ndim, rect, niter, true);
    norm = 0.;
    for (id=0; id < nd; id++) {
	a = cabsf(c_den[id]);
	norm += a*a;
    }
    norm = sqrtf(nd/norm);

    for (id=0; id < nd; id++) {
#ifdef SF_HAS_COMPLEX_H
	c_num[id] *= norm;
	c_den[id] *= norm;
#else
	c_num[id] = sf_crmul(c_num[id],norm);
	c_den[id] = sf_crmul(c_den[id],norm);
#endif
    }
    /* output the ratio to w_freq */
    cdivn (c_num, c_den, cw_freq);
    /* Convert to real */
    for(id=0; id<nd;id++){
	w_freq[id] =crealf(cw_freq[id]);
    }

}
