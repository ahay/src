#include <math.h>

#include <rsf.h>

#include "monofshape.h"
#include "monof.h"

static int nfft, nw;
static float complex *cdata;
static float *shape, *tmp, dw;
static kiss_fftr_cfg forw, invs;

void monofshape_init(int n1)
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = n1;
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    cdata = sf_complexalloc(nw);
    shape = sf_floatalloc(nw);
    tmp = sf_floatalloc(nfft);
    
    forw = kiss_fftr_alloc(nfft,0,NULL,NULL);
    invs = kiss_fftr_alloc(nfft,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation problem",__FILE__);
}

void monofshape_set(float a0, int n, float* pattern, int niter)
{
    int iw, i0;
    float f, w, max, scale;
    
    for (iw=0; iw < n; iw++) {
	tmp[iw] = pattern[iw];
    }
    for (iw=n; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    scale = sqrtf(1./nfft); /* FFT scaling */ 

    kiss_fftr(forw, tmp, (kiss_fft_cpx *) cdata);
    max = 0.;
    i0 = 0;
    for (iw=0; iw < nw; iw++) {
	f = cabsf(cdata[iw])*scale;
	if (f > max) {
	    max = f;
	    i0 = iw;
	}
	tmp[iw] = f;
    }

    sf_warning("i0=%d max=%g",i0,max);

    f = monof(tmp,i0,niter,a0,nw,dw,true);

    for (iw=0; iw < nw; iw++) {
	w = (iw-i0)*dw;
	w *= w;
	shape[iw] = expf(-0.5*f*w)/nfft;
    } 
    
}

void monofshape_close(void) {
    free(cdata);
    free(shape);
    free(tmp);
    free(forw);
    free(invs);
}

void monofshape_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
{
    int iw;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (iw=0; iw < nx; iw++) {	    
	tmp[iw] = adj? y[iw] : x[iw];
    }

    for (iw = nx; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    kiss_fftr(forw, tmp, (kiss_fft_cpx *) cdata);
    for (iw=0; iw < nw; iw++) {
	cdata[iw] *= shape[iw];
    }
    kiss_fftri(invs,(const kiss_fft_cpx *) cdata, tmp);

    for (iw = 0; iw < nx; iw++) {
	if (adj) {
	    x[iw] += tmp[iw];
	} else {
	    y[iw] += tmp[iw];
	}
    }
}

/* 	$Id$	 */
