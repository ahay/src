#include <math.h>

#include <rsf.h>

#include "monofshape.h"
#include "monof.h"

static int nfft, nw;
static float complex *cdata;
static float *shape, *tmp, dw;

void monofshape_init(int n1)
{
    /* determine frequency sampling (for real to complex FFT) */
    nfft = sf_npfar(n1);
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    cdata = sf_complexalloc(nw);
    shape = sf_floatalloc(nw);
    tmp = sf_floatalloc(nfft);
}

void monofshape_set(float a0, int n, float* pattern, int niter)
{
    int iw, i0;
    float f, w, max;
    
    for (iw=0; iw < n; iw++) {
	tmp[iw] = pattern[iw];
    }
    for (iw=n; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    sf_pfarc (1,nfft,tmp,cdata);
    max = 0.;
    i0 = 0;
    for (iw=0; iw < nw; iw++) {
	f = cabsf(cdata[iw]);
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
	shape[iw] = expf(-f*w)/nfft;
    } 
}

void monofshape_close(void) {
    free(cdata);
    free(shape);
    free(tmp);
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

    sf_pfarc (1,nfft,tmp,cdata);
    for (iw=0; iw < nw; iw++) {
	cdata[iw] *= shape[iw];
    }
    sf_pfacr (-1,nfft,cdata,tmp);	

    for (iw = 0; iw < nx; iw++) {
	if (adj) {
	    x[iw] += tmp[iw];
	} else {
	    y[iw] += tmp[iw];
	}
    }
}

/* 	$Id: monofshape.c,v 1.1 2004/02/14 06:57:16 fomels Exp $	 */
