#include <math.h>

#include <rsf.h>

#include "freqfilt2.h"

static int nfft, nw, nk, m1, m2;
static float complex **cdata, *ctrace;
static float **shape, **tmp;

void freqfilt2_init(int n1, int n2, int nfft1, int nw1, int nk1)
{
    m1 = n1;
    nw = nw1;
    m2 = n2;
    nfft = nfft1;
    nk = nk1;

    cdata = sf_complexalloc2(n1,nk);
    ctrace = sf_complexalloc(nw);
    tmp = sf_floatalloc2(n1,nfft);
}

void freqfilt2_set(float **filt)
{
    shape = filt;
}

void freqfilt2_close(void) {
    free(cdata[0]);
    free(cdata);
    free(tmp[0]);
    free(tmp);
    free(ctrace);
}

void freqfilt2_spec (const float* x, float** y) {
    int ik, iw;

    for (ik=0; ik < m2; ik++) {
	for (iw=0; iw < m1; iw++) {
	    tmp[ik][iw] = x[ik*m1+iw];
	}
    }

    for (ik=m2; ik < nfft; ik++) {
	for (iw=0; iw < m1; iw++) {
	    tmp[ik][iw] = 0.;
	}
    }
    	
    sf_pfa2rc (1,2,m1,nfft,tmp[0],cdata[0]);

    for (ik=0; ik < nk; ik++) {
	for (iw=0; iw < m1; iw++) {
	    ctrace[iw] = iw%2? -cdata[ik][iw]: cdata[ik][iw];
	}
	for (iw=m1; iw < nw; iw++) {
	    ctrace[iw] = 0.;
	}
	
	sf_pfacc (1,nw,ctrace);

	for (iw=0; iw < nw; iw++) {
	    y[ik][iw] = cabsf(ctrace[iw]);
	}
    }
}

void freqfilt2_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
{
    int iw, ik;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (ik=0; ik < m2; ik++) {
	for (iw=0; iw < m1; iw++) {
	    tmp[ik][iw] = adj? y[ik*m1+iw]: x[ik*m1+iw];
	}
    }

    for (ik=m2; ik < nfft; ik++) {
	for (iw=0; iw < m1; iw++) {
	    tmp[ik][iw] = 0.;
	}
    }
    	
    sf_pfa2rc (1,2,m1,nfft,tmp[0],cdata[0]);

    for (ik=0; ik < nk; ik++) {
	for (iw=0; iw < m1; iw++) {
	    ctrace[iw] = iw%2? -cdata[ik][iw]: cdata[ik][iw];
	}
	for (iw=m1; iw < nw; iw++) {
	    ctrace[iw] = 0.;
	}

	sf_pfacc (1,nw,ctrace);

	for (iw=0; iw < nw; iw++) {
	    ctrace[iw] *= shape[ik][iw];
	}

	sf_pfacc (-1,nw,ctrace);

	for (iw=0; iw < m1; iw++) {
	    cdata[ik][iw] = iw%2? -ctrace[iw]: ctrace[iw];
	}
    }

    sf_pfa2cr (-1,2,m1,nfft,cdata[0],tmp[0]);

    for (ik=0; ik < m2; ik++) {
	for (iw=0; iw < m1; iw++) {	  
	    if (adj) {
		x[ik*m1+iw] += tmp[ik][iw];
	    } else {
		y[ik*m1+iw] += tmp[ik][iw];
	    }
	}
    }
}

/* 	$Id: freqfilt2.c,v 1.1 2004/02/26 05:20:50 fomels Exp $	 */
