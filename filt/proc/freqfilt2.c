#include <math.h>

#include <rsf.h>

#include "freqfilt2.h"

static int nfft, nw, m1, m2;
static float complex *ctrace, *ctrace2, **fft;
static float *trace, **shape;
kiss_fftr_cfg tfor, tinv;
kiss_fft_cfg  xfor, xinv;

void freqfilt2_init(int n1, int n2, int nw1)
{
    m1 = n1;
    nw = nw1;
    m2 = n2;
    nfft = n1;
    if (n1%2) nfft++;

    tfor = kiss_fftr_alloc(nfft,0,NULL,NULL);
    tinv = kiss_fftr_alloc(nfft,1,NULL,NULL);
    xfor = kiss_fft_alloc(n2,0,NULL,NULL);
    xinv = kiss_fft_alloc(n2,1,NULL,NULL);
    if (NULL == tfor || NULL == tinv || NULL == xfor || NULL == xinv)
	sf_error("%s: KISS FFT allocation error",__FILE__);

    trace = sf_floatalloc(nfft);
    ctrace = sf_complexalloc(nw);
    ctrace2 = sf_complexalloc(n2);
    fft = sf_complexalloc2(nw,n2);

    if (n1%2) trace[n1]=0.;
}

void freqfilt2_set(float **filt)
{
    shape = filt;
}

void freqfilt2_close(void) {
    free (tfor);
    free (tinv);
    free (xfor);
    free (xinv);
    free (trace);
    free (ctrace);
    free (ctrace2);
    free (fft[0]);
    free (fft);
}

void freqfilt2_spec (const float* x, float** y) {
    int ik, iw;

    for (ik=0; ik < m2; ik++) {
	for (iw=0; iw < m1; iw++) {
	    trace[iw] = x[ik*m1+iw];
	}
	kiss_fftr (tfor,trace,(kiss_fft_cpx *) ctrace);
	for (iw=0; iw < nw; iw++) {
	    fft[ik][iw] = ik%2? -ctrace[iw]: ctrace[iw];
	}
    }

    for (iw=0; iw < nw; iw++) {
	kiss_fft_stride(xfor,(kiss_fft_cpx *) (fft[0]+iw),
			(kiss_fft_cpx *) ctrace2,nw);
	for (ik=0; ik < m2; ik++) {
	    y[iw][ik] = cabsf(ctrace2[ik]); /* transpose */
	}
    }
}

void freqfilt2_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
{
    int iw, ik;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (ik=0; ik < m2; ik++) {
	for (iw=0; iw < m1; iw++) {
	    trace[iw] = adj? y[ik*m1+iw]: x[ik*m1+iw];
	}
	kiss_fftr (tfor,trace,(kiss_fft_cpx *) ctrace);
	for (iw=0; iw < nw; iw++) {
	    fft[ik][iw] = ik%2? -ctrace[iw]: ctrace[iw];
	}
    }

    for (iw=0; iw < nw; iw++) {
	kiss_fft_stride(xfor,(kiss_fft_cpx *) (fft[0]+iw),
			(kiss_fft_cpx *) ctrace2,nw);

	for (ik=0; ik < m2; ik++) {
	    ctrace2[ik] *= shape[iw][ik];
	}

	kiss_fft(xinv,(kiss_fft_cpx *) ctrace2,(kiss_fft_cpx *) ctrace2);

	for (ik=0; ik < m2; ik++) {
	    fft[ik][iw] = ik%2? -ctrace2[ik]: ctrace2[ik];
	}
    }

    for (ik=0; ik < m2; ik++) {
	kiss_fftri (tinv,(kiss_fft_cpx *) fft[ik], trace);

	for (iw=0; iw < m1; iw++) {	  
	    if (adj) {
		x[ik*m1+iw] += trace[iw];
	    } else {
		y[ik*m1+iw] += trace[iw];
	    }
	}
    }
}

/* 	$Id$	 */
