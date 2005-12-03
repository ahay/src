/* Frequency spectra in 2-D.
*/
/*
Copyright (C) 2004 University of Texas at Austin

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[]) 
{
    int nw, n1, n2, n3, ni, nfft, i, i1, i2, i3;
    float d1, o1, d2, o2, dw, dk, k0, **spec, scale, *trace;
    float complex **fft, *ctrace, *ctrace2;
    char key[3];
    bool sum;
    kiss_fftr_cfg tfft;
    kiss_fft_cfg  xfft;
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
    nfft = sf_fftr_size(n1);
    if (nfft%2) nfft++;
    nw = nfft/2+1;
    dw = 1./(nfft*d1);

    /* determine wavenumber sampling (for real to complex FFT) */
    dk = 1./(n2*d2);
    k0 = -0.5/d2;

    trace = sf_floatalloc (nfft);
    ctrace = sf_complexalloc (nw);
    ctrace2 = sf_complexalloc (n2);
    fft = sf_complexalloc2(nw,n2);
    spec = sf_floatalloc2(nw,n2);
 	
    tfft = kiss_fftr_alloc(nfft,0,NULL,NULL);
    xfft = kiss_fft_alloc(n2,0,NULL,NULL);

    sf_putint(out,"n1",nw);
    sf_putfloat(out,"d1",dw);
    sf_putfloat(out,"o1",0.);

    sf_putfloat(out,"d2",dk);
    sf_putfloat(out,"o2",k0);

    if (sum) {
	for (i=2; i < SF_MAX_DIM; i++) {
	    snprintf(key,3,"n%d",i+1);
	    if (!sf_histint(in,key,&ni)) break;
	    if (ni > 1) sf_putint(out,key,1);
	}
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < nw; i1++) {
		spec[i2][i1] = 0.;
	    }
	}
    }
    
    scale = sqrtf(1./(nfft*n2)); /* FFT scaling */ 

    for (i1=n1; i1 < nfft; i1++) {
	trace[i1]=0.;
    }

    /*  loop over all planes */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread(trace,n1,in);
	    /* Fourier transform in n1 */
	    kiss_fftr (tfft,trace,(kiss_fft_cpx *) ctrace);
	    for (i1=0; i1 < nw; i1++) {
		fft[i2][i1] = i2%2? -ctrace[i1]: ctrace[i1];
	    }
	}

	for (i1=0; i1 < nw; i1++) {
	    /* Fourier transform in n2 */
	    kiss_fft_stride(xfft,(kiss_fft_cpx *) (fft[0]+i1),
			    (kiss_fft_cpx *) ctrace2,nw);
	    if (sum) {
		for (i2=0; i2 < n2; i2++) {
		    spec[i2][i1] += cabsf(ctrace2[i2]);
		}
	    } else {
		for (i2=0; i2 < n2; i2++) {
		    spec[i2][i1] = cabsf(ctrace2[i2])*scale;
		}
	    }
	}
	
	if (!sum) sf_floatwrite(spec[0],nw*n2,out);
    } /* i3 */

    if (sum) { /* normalize and output */
	scale /= n3;
	for (i2=0; i2 < n2; i2++) { 
	    for (i1=0; i1 < nw; i1++) {
		spec[i2][i1] *= scale;
	    }
	}
	sf_floatwrite(spec[0],nw*n2,out);
    }

    exit (0);
}

/* 	$Id$	 */
