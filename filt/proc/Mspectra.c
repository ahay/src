/* Frequency spectra.

Takes: < data.rsf > spectra.rsf
*/

#include <stdio.h>
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[]) 
{
    int n1, n2, ni, nfft, nw, i, i1, i2;
    float d1, o1, dw, *spec, *trace, scale;
    complex float *fft;
    char key[3];
    bool sum, isphase;
    sf_file in, out;

    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype (in)) sf_error("Need float data");

    if (!sf_histint(in,"n1",&n1)) n1=1;
    n2 = sf_leftsize(in,1);

    if (!sf_getbool("all",&sum)) sum=false;
    /* if y, compute average spectrum for all traces */

    if (!sf_histfloat(in,"d1",&d1)) d1=0.004;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    /* determine frequency sampling (for real to complex FFT) */
    nfft = sf_npfar(n1);
    nw = nfft/2+1;
    dw = 1./(nfft*d1);

    trace = sf_floatalloc (nfft);
    fft = sf_complexalloc (nw);
    spec = sf_floatalloc (nw);
 	
    sf_putint(out,"n1",nw);
    sf_putfloat(out,"d1",dw);
    sf_putfloat(out,"o1",0.);

    if (sum) {
	for (i=1; i < SF_MAX_DIM; i++) {
	    snprintf(key,3,"n%d",i+1);
	    if (!sf_histint(in,key,&ni)) break;
	    if (ni > 1) sf_putint(out,key,1);
	}
	for (i1=0; i1 < nw; i1++) {
	    spec[i1] = 0.;
	}
    }

    if (!sf_getbool("phase",&isphase)) isphase=false;
    /* if y, compute phase spectra */
    
    for (i1=n1; i1 < nfft; i1++) { /* pad with zeros */
	trace[i1]=0.;
    }

    scale = sqrtf(1./nfft); /* FFT scaling */ 

/*  loop over all traces */
    for (i2=0; i2 < n2; i2++) {
	sf_read(trace,sizeof(float),n1,in);

	/* Fourier transform */
	sf_pfarc (1,nfft,trace,fft);

	if (sum) {
	    if (isphase) {
		for (i1=0; i1 < nw; i1++) {
		    spec[i1] += cargf(fft[i1]*cexpf(2.*I*SF_PI*i1*dw*o1));
		}
	    } else {
		for (i1=0; i1 < nw; i1++) {
		    spec[i1] += cabsf(fft[i1]);
		}
	    }
	} else {
	    if (isphase) {
		for (i1=0; i1 < nw; i1++) {
		    spec[i1] = cargf(fft[i1]*cexpf(2.*I*SF_PI*i1*dw*o1));
		}
	    } else {
		for (i1=0; i1 < nw; i1++) {
		    spec[i1] = cabsf(fft[i1])*scale;
		}
	    }
	    sf_write(spec,sizeof(float),nw,out);
	}
    }

    if (sum) { /* normalize and output */
	scale /= n2;
	for (i1=0; i1 < nw; i1++) {
	    spec[i1] *= scale;
	}
	sf_write(spec,sizeof(float),nw,out);
    }

    exit (0);
}

/* 	$Id: Mspectra.c,v 1.2 2003/10/01 22:45:55 fomels Exp $	 */

