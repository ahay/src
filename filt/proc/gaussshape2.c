#include <math.h>

#include <rsf.h>

#include "monof2.h"
#include "gaussshape2.h"

static int nfft, nw, nk, m1, m2;
static float complex **cdata, *ctrace;
static float **shape, **tmp, dw, dk, w0;

void gaussshape2_init(int n1, int n2)
{
    /* determine frequency sampling (for real to complex FFT) */
    m1 = n1;
    nw = sf_npfa(n1);
    dw = 2.*SF_PI/nw;
    w0 = -SF_PI;

    /* determine wavenumber sampling (for complex FFT) */
    m2 = n2;
    nfft = sf_npfar(n2);
    nk = nfft/2+1;
    dk = 2.*SF_PI/nfft;

    cdata = sf_complexalloc2(n1,nk);
    ctrace = sf_complexalloc(nw);
    shape = sf_floatalloc2(nw,nk);
    tmp = sf_floatalloc2(n1,nfft);
}

void gaussshape2_set(float* a, int n1, int n2, float* pattern, int niter)
{
    int ik, iw;
    float w, w2, k, k2, wk;
    
    for (ik=0; ik < n2; ik++) {
	for (iw=0; iw < n1; iw++) {
	    tmp[ik][iw] = pattern[ik*n1+iw];
	}
    }

    for (ik=n2; ik < nfft; ik++) {
	for (iw=0; iw < n1; iw++) {
	    tmp[ik][iw] = 0.;
	}
    }
    	
    sf_pfa2rc (1,2,n1,nfft,tmp[0],cdata[0]);

    for (ik=0; ik < nk; ik++) {
	for (iw=0; iw < n1; iw++) {
	    ctrace[iw] = iw%2? -cdata[ik][iw]: cdata[ik][iw];
	}
	for (iw=n1; iw < nw; iw++) {
	    ctrace[iw] = 0.;
	}
	
	sf_pfacc (1,nw,ctrace);
	
	for (iw=0; iw < nw; iw++) {
	    shape[ik][iw] = cabsf(ctrace[iw]);
	}
    }

    monof2(shape, niter, a, nw, dw, w0, nk, dk, 0., true);

    for (ik=0; ik < nk; ik++) {
	k = ik*dk;
	k2 = k*k;
	for (iw=0; iw < nw; iw++) {
	    w = w0+iw*dw;
	    w2 = w*w;
	    wk = w*k;
	    shape[ik][iw] = expf(-a[0]*w2-a[1]*wk-a[2]*k2)/(nfft*nw);
	}
    }
}

void gaussshape2_close(void) {
    free(cdata[0]);
    free(cdata);
    free(shape[0]);
    free(shape);
    free(tmp[0]);
    free(tmp);
    free(ctrace);
}

void gaussshape2_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
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

    sf_pfa2cr (1,2,m1,nfft,cdata[0],tmp[0]);

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

/* 	$Id: gaussshape2.c,v 1.1 2004/02/25 16:16:27 fomels Exp $	 */
