#include <math.h>

#include <rsf.h>

#include "gauss.h"

static int nfft, nw;
static float complex *cdata;
static float *shape, *tmp;

void gauss_init(int n1, float rect)
{
    int iw;
    float dw, w;

    /* determine frequency sampling (for real to complex FFT) */
    nfft = sf_npfar(n1);
    nw = nfft/2+1;
    dw = 2.*SF_PI/nfft;

    cdata = sf_complexalloc(nw);
    shape = sf_floatalloc(nw);
    tmp = sf_floatalloc(nfft);

    rect = sqrtf((rect*rect-1.)/12.);

    for (iw=0; iw < nw; iw++) {
	w = iw*dw*rect;
	w *= w;
	shape[iw] = expf(-w)/nfft;
    }
}

void gauss_close(void) {
    free(cdata);
    free(shape);
    free(tmp);
}

void gauss_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
{
    int iw;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (iw=0; iw < nx; iw++) {
	tmp[iw] = adj? y[iw] : x[iw];
    }
    for (iw=nx; iw < nfft; iw++) {
	tmp[iw] = 0.;
    }

    sf_pfarc (1,nfft,tmp,cdata);

    for (iw=0; iw < nw; iw++) {
	cdata[iw] *= shape[iw];
    }
    sf_pfacr (-1,nfft,cdata,tmp);	

    for (iw=0; iw < nx; iw++) {	    
	if (adj) {
	    x[iw] += tmp[iw];
	} else {
	    y[iw] += tmp[iw];
	}
    } 
}

/* 	$Id: gauss.c,v 1.1 2004/02/14 06:59:24 fomels Exp $	 */
