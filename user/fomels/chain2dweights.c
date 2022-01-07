/* space and frequency weights of 2D chain */
#include <rsf.h>
#include "chain2dweights.h"
#include "fft2.h"


static int nz, nx, nz1, nx2, nk_chain;
static float *w, *wf, **tmp1;
static sf_complex *ctmp1; /* for 2D-fft */

void chain2dweights_init(int n1, /* trace length */
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
    nk_chain = nk_in;
    nz1 = nz2_in;
    nx2 = nx2_in;
    w = ww;
    wf = ff;
    /*for 2D fft*/
    tmp1 = sf_floatalloc2(nz1,nx2);  
    ctmp1 = sf_complexalloc(nk_chain);
}

void chain2dweights_close(void)
/*< clean allocated storage >*/
{
    free(*tmp1);
    free(tmp1);
    free(ctmp1);
}


void chain2dweights_lop(bool adj, bool add, int nxx, int nyy, float* x, float* y) 
/*< linear operator >*/
{


	int i, i1, i2, ik;
    if (nxx != nz*nx || nyy != nz*nx) sf_error("%s: Wrong size",__FILE__);


    sf_adjnull(adj,add,nxx,nyy,x,y);

    if (adj){ /* Adjoint*/

	   /* pad with zeros */
	    for (i2=0; i2 < nx; i2++) {
		for (i1=0; i1 < nz; i1++) {
		    i = i1+i2*nz;
		    tmp1[i2][i1] = w[i]*y[i];
		}
		for (i1=nz; i1 < nz1; i1++) {
		    tmp1[i2][i1] = 0.0f;
		}
	    }
	    for (i2=nx; i2 < nx2; i2++) {
		for (i1=0; i1 < nz1; i1++) {
		    tmp1[i2][i1] = 0.0f;
		}
	    }

	    /* forward FFT */
	    fft2_allocate(ctmp1);
	    fft2(tmp1[0],ctmp1);

	    /* frequency weight */
	    for (ik=0; ik < nk_chain; ik++) {
		ctmp1[ik] *= wf[ik];
	    }
	    /* inverse FFT */
	    ifft2(tmp1[0],ctmp1);

	    for (i=0; i < nz*nx; i++) {
		i1 = i%nz;
		i2 = i/nz;
		x[i] += tmp1[i2][i1];
	    }

    } else{ /* Forward */

	   /* pad with zeros */
	    for (i2=0; i2 < nx; i2++) {
		for (i1=0; i1 < nz; i1++) {
		    i = i1+i2*nz;
		    tmp1[i2][i1] = x[i];
		}
		for (i1=nz; i1 < nz1; i1++) {
		    tmp1[i2][i1] = 0.0f;
		}
	    }
	    for (i2=nx; i2 < nx2; i2++) {
		for (i1=0; i1 < nz1; i1++) {
		    tmp1[i2][i1] = 0.0f;
		}
	    }

	    /* forward FFT */
	    fft2_allocate(ctmp1);
	    fft2(tmp1[0],ctmp1);

	    /* frequency weight */
	    for (ik=0; ik < nk_chain; ik++) {
		ctmp1[ik] *= wf[ik];
	    }
	    /* inverse FFT */
	    ifft2(tmp1[0],ctmp1);

	    for (i=0; i < nz*nx; i++) {
		i1 = i%nz;
		i2 = i/nz;
		y[i] += w[i]*tmp1[i2][i1];
	    }
    }
}
