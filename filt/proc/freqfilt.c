#include <math.h>

#include <rsf.h>

#include "freqfilt.h"

static int nfft, nw, n;
static float complex *cdata;
static float *shape, *tmp;

void freqfilt_init(int n1, int nfft1, int nw1)
{
    n = n1;
    nfft = nfft1;
    nw = nw1;

    cdata = sf_complexalloc(nw);
    shape = sf_floatalloc(nw);
    tmp = sf_floatalloc(nfft);
}

void freqfilt_set(float *filt)
{
    shape = filt;
}

void freqfilt_close(void) {
    free(cdata);
    free(tmp);
}

void freqfilt_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
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

/* 	$Id: freqfilt.c,v 1.1 2004/04/02 02:30:49 fomels Exp $	 */
