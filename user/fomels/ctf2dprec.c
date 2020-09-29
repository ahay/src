/* TF Weights Preconditioner for Complex input as linear oper. */
#include <rsf.h>
#include "ctf2dprec.h"
#include "cfft2w.h"


static int nz, nx, nz2, nx2, nk;
static float *w, *wf;
static sf_complex **tmp1, *ctmp1; /* for 2D-fft */



void ctf2dprec_init(int n1, /* trace length */
			int n2, /* number of length */
			int nk_in, /* total Fourier size [n1*n2] */
			int nz2_in, /* Fourier dim1 */
			int nx2_in, /* Fourier dim2 */
		   float* ww /* [n1*n2] time weight */,
		   float* ff /* [nk] frequency weight */)
/*< initialize >*/

{
    nz = n1;
    nx = n2;
	nk = nk_in;
    nz2 = nz2_in;
    nx2 = nx2_in;
    w = ww;
    wf = ff;
    /*for 2D fft*/
    tmp1 = sf_complexalloc2(nz2,nx2);  
    ctmp1 = sf_complexalloc(nk);
}


void ctf2dprec_lop(bool adj, bool add, int nxx, int nyy, sf_complex* x, sf_complex* y) 
/*< linear operator >*/
{

	int i, i1, i2, ik;
    if (nxx != nz*nx || nyy != nz*nx) sf_error("%s: Wrong size",__FILE__);

    sf_cadjnull(adj,add,nxx,nyy,x,y);

    if (adj){ /* Adjoint*/

	   /* pad with zeros */
	    for (i2=0; i2 < nx; i2++) {
		for (i1=0; i1 < nz; i1++) {
		    i = i1+i2*nz;
	#ifdef SF_HAS_COMPLEX_H
		    tmp1[i2][i1] = w[i]*y[i];
	#else
		    tmp1[i2][i1] = sf_crmul(w[i],y[i]);
	#endif
		}
		for (i1=nz; i1 < nz2; i1++) {
		    tmp1[i2][i1] = sf_cmplx(0.0,0.0);
		}
	    }
	    for (i2=nx; i2 < nx2; i2++) {
		for (i1=0; i1 < nz2; i1++) {
		    tmp1[i2][i1] = sf_cmplx(0.0,0.0);
		}
	    }

	    /* forward FFT */
	    icfft2_allocate(ctmp1);
	    cfft2(tmp1[0],ctmp1);

	    /* frequency weight */
/*	    for (ik=0; ik < nk; ik++) {
	#ifdef SF_HAS_COMPLEX_H
			ctmp1[ik] *= wf[ik];
	#else
		    ctmp1[ik] = sf_crmul(wf[ik],ctmp1[ik]);
	#endif
	    }*/
	    /* inverse FFT */
	    icfft2(tmp1[0],ctmp1);

	    for (i=0; i < nz*nx; i++) {
		i1 = i%nz;
		i2 = i/nz;
	#ifdef SF_HAS_COMPLEX_H
		x[i] += tmp1[i2][i1];
	#else
		x[i] = sf_cadd(x[i],tmp1[i2][i1]);
	#endif
	    }
	cfft2_finalize();

    } else{ /* Forward */

	   /* pad with zeros */
	    for (i2=0; i2 < nx; i2++) {
		for (i1=0; i1 < nz; i1++) {
		    i = i1+i2*nz;
		    tmp1[i2][i1] = x[i];
		}
		for (i1=nz; i1 < nz2; i1++) {
		    tmp1[i2][i1] = sf_cmplx(0.0,0.0);
		}
	    }
	    for (i2=nx; i2 < nx2; i2++) {
		for (i1=0; i1 < nz2; i1++) {
		    tmp1[i2][i1] = sf_cmplx(0.0,0.0);
		}
	    }
	    /* forward FFT */
	    icfft2_allocate(ctmp1);
	    cfft2(tmp1[0],ctmp1);

	    /* frequency weight */
/*	    for (ik=0; ik < nk; ik++) {
	#ifdef SF_HAS_COMPLEX_H
			ctmp1[ik] *= wf[ik];
	#else
		    ctmp1[ik] = sf_crmul(wf[ik],ctmp1[ik]);
	#endif
	    }*/
	    /* inverse FFT */
	    icfft2(tmp1[0],ctmp1);

	    for (i=0; i < nz*nx; i++) {
		i1 = i%nz;
		i2 = i/nz;
	#ifdef SF_HAS_COMPLEX_H
		y[i] += w[i]*tmp1[i2][i1];
	#else
		y[i] = sf_cadd(y[i],sf_crmul(w[i],tmp1[i2][i1]));
	#endif
	    }

		cfft2_finalize();
    }
}
