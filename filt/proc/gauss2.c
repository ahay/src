#include <math.h>

#include <rsf.h>

#include "gauss2.h"

static int nfft, nw, nk, m1, m2;
static float complex **cdata, *ctrace;
static float **shape, **tmp;

void gauss2_init(int n1, int n2, float f1, float f2)
{
    int ik, iw;
    float dw, w, w2, dk, k, w0, k2;

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

    f1 = sqrtf((f1*f1-1.)/12.);
    f2 = sqrtf((f2*f2-1.)/12.);

    for (ik=0; ik < nk; ik++) {
	k = ik*dk;
	k2 = f2*k;
	k2 *= k2;
	for (iw=0; iw < nw; iw++) {
	    w = w0+iw*dw;
	    w2 = f1*w;
	    w2 *= w2;
	    shape[ik][iw] = expf(-k2-w2)/(nfft*nw);
	}
    }
}

void gauss2_close(void) {
    free(cdata[0]);
    free(cdata);
    free(shape[0]);
    free(shape);
    free(tmp[0]);
    free(tmp);
    free(ctrace);
}

void gauss2_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
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

/* 	$Id: gauss2.c,v 1.1 2004/02/14 06:57:16 fomels Exp $	 */
