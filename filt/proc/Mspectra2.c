/* Frequency spectra in 2-D.

Takes: < data.rsf > spectra.rsf
*/

#include <stdio.h>
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[]) 
{
    int n1, n2, n3, ni, nfft, nw, nk, i, i1, i2, i3;
    float d1, o1, d2, o2, dw, dk, w0, **spec, scale, **plane;
    float complex **fft, *trace;
    char key[3];
    bool sum;
    sf_file in, out;

    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype (in)) sf_error("Need float data");

    if (!sf_histint(in,"n1",&n1)) n1=1;
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    if (!sf_getbool("all",&sum)) sum=false;
    /* if y, compute average spectrum for all traces */

    if (!sf_histfloat(in,"d1",&d1)) d1=0.004;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    /* determine frequency sampling (for complex FFT) */
    nw = sf_npfa(n1);
    dw = 1./(nw*d1);
    w0 = -0.5/d1;

    /* determine wavenumber sampling (for real to complex FFT) */
    nfft = sf_npfar(n2);
    nk = nfft/2+1;
    dk = 1./(nfft*d2);

    trace = sf_complexalloc (nw);
    plane = sf_floatalloc2 (n1,nfft);
    fft = sf_complexalloc2 (n1,nk);
    spec = sf_floatalloc2 (nw,nk);
 	
    sf_putint(out,"n1",nw);
    sf_putfloat(out,"d1",dw);
    sf_putfloat(out,"o1",w0);

    sf_putint(out,"n2",nk);
    sf_putfloat(out,"d2",dk);
    sf_putfloat(out,"o2",0.);

    if (sum) {
	for (i=2; i < SF_MAX_DIM; i++) {
	    snprintf(key,3,"n%d",i+1);
	    if (!sf_histint(in,key,&ni)) break;
	    if (ni > 1) sf_putint(out,key,1);
	}
	for (i2=0; i2 < nk; i2++) {
	    for (i1=0; i1 < nw; i1++) {
		spec[i2][i1] = 0.;
	    }
	}
    }
    
    for (i2=n2; i2 < nfft; i2++) { /* pad with zeros */
	for (i1=0; i1 < n1; i1++) {
	    plane[i2][i1]=0.;
	}
    }

    scale = sqrtf(1./(nfft*nw)); /* FFT scaling */ 

    /*  loop over all planes */
    for (i3=0; i3 < n3; i3++) {
	sf_read(plane[0],sizeof(float),n1*n2,in);

	/* Fourier transform */
	sf_pfa2rc (1,2,n1,nfft,plane[0],fft[0]);

	for (i2=0; i2 < nk; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = i1%2? -fft[i2][i1]: fft[i2][i1];
	    }
	    for (i1=n1; i1 < nw; i1++) {
		trace[i1] = 0.;
	    }

	    sf_pfacc (1,nw,trace);

	    if (sum) {
		for (i1=0; i1 < nw; i1++) {
		    spec[i2][i1] += cabsf(trace[i1]);
		}
	    } else {
		for (i1=0; i1 < nw; i1++) {
		    spec[i2][i1] = cabsf(trace[i1])*scale;
		}
	    }
	} /* i2 */

	sf_write(spec[0],sizeof(float),nw*nk,out);
    } /* i3 */

    if (sum) { /* normalize and output */
	scale /= n3;
	for (i2=0; i2 < nk; i2++) { 
	    for (i1=0; i1 < nw; i1++) {
		spec[i2][i1] *= scale;
	    }
	}
	sf_write(spec[0],sizeof(float),nw*nk,out);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mspectra2.c,v 1.2 2004/03/22 05:43:25 fomels Exp $	 */

