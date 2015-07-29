/* Frequency spectra in 2-D. */
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
#include <rsf.h>

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

int main (int argc, char* argv[]) 
{
    int nw, n1, n2, nk, n3, ni, nfft, i, i1, i2, i3;
    float d1, o1, d2, o2, dw, dk, k0, **spec, scale, *trace, **traces;
    kiss_fft_cpx **fft;
    char key[3], *label;
    bool sum;
    sf_file in, out;

#ifdef SF_HAS_FFTW
    fftwf_plan txfft;
#else
    kiss_fft_cpx *ctrace, *ctrace2;
    kiss_fftr_cfg tfft;
    kiss_fft_cfg  xfft;
#endif

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

    /* determine frequency sampling (for real to complex FFT) */
    nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
    nw = nfft/2+1;
    dw = 1./(nfft*d1);

    /* determine wavenumber sampling (for complex FFT) */
    nk = kiss_fft_next_fast_size(n2);
    dk = 1./(nk*d2);
    k0 = -0.5/d2;

    traces = sf_floatalloc2 (nfft,nk);
    fft = (kiss_fft_cpx**) sf_complexalloc2(nw,nk);
    spec = sf_floatalloc2(nw,nk);

#ifdef SF_HAS_FFTW
    txfft = fftwf_plan_dft_r2c_2d(nk,nfft,
				  traces[0], (fftwf_complex *) fft[0],
				  FFTW_ESTIMATE);
    if (NULL == txfft) sf_error("FFTW failure.");
#else
    ctrace = (kiss_fft_cpx*) sf_complexalloc (nw);
    ctrace2 = (kiss_fft_cpx*) sf_complexalloc (nk);

    tfft = kiss_fftr_alloc(nfft,0,NULL,NULL);
    xfft = kiss_fft_alloc(nk,0,NULL,NULL);
#endif

    sf_putint(out,"n1",nw);
    sf_putfloat(out,"d1",dw);
    sf_putfloat(out,"o1",0.);

    sf_putint(out,"n2",nk);
    sf_putfloat(out,"d2",dk);
    sf_putfloat(out,"o2",k0);

    /* fix labels and units */
    if (NULL != (label = sf_histstring(in,"label1")) &&
	!sf_fft_label(1,label,out))
	sf_putstring(out,"label1","Wavenumber");
    if (NULL != (label = sf_histstring(in,"label2")) &&
	!sf_fft_label(2,label,out))
	sf_putstring(out,"label2","Wavenumber");

    sf_fft_unit(1,sf_histstring(in,"unit1"),out);
    sf_fft_unit(2,sf_histstring(in,"unit2"),out);

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

    scale = sqrtf(1./(nfft*nk)); /* FFT scaling */

    for (i2=n2; i2 < nk; i2++) {
	for (i1=0; i1 < nfft; i1++) {
	    traces[i2][i1] = 0.;
	}
    }

    /*  loop over all planes */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    trace = traces[i2];
	    
	    sf_floatread(trace,n1,in);

	    if (i2%2) { /* FFT centering */
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] = -trace[i1];
		}
	    }
	    for (i1=n1; i1 < nfft; i1++) {
		trace[i1]=0.;
	    }

#ifndef SF_HAS_FFTW
	    /* Fourier transform in n1 */
	    kiss_fftr (tfft,trace,ctrace);
	    for (i1=0; i1 < nw; i1++) {
		fft[i2][i1] = ctrace[i1];
	    }
#endif
	}
#ifdef SF_HAS_FFTW

	fftwf_execute(txfft);

	for (i2=0; i2 < nk; i2++) {
	    for (i1=0; i1 < nw; i1++) {
		if (sum) {
		    spec[i2][i1] += sf_cabsf(fft[i2][i1]);
		} else { 
		    spec[i2][i1] = sf_cabsf(fft[i2][i1])*scale;
		}
	    }
	}
#else
	for (i2=n2; i2 < nk; i2++) {
	    for (i1=0; i1 < nw; i1++) {
		fft[i2][i1].r = 0.;
		fft[i2][i1].i = 0.;
	    }
	}

	for (i1=0; i1 < nw; i1++) {
	    /* Fourier transform in n2 */
	    kiss_fft_stride(xfft,fft[0]+i1,ctrace2,nw);
	    for (i2=0; i2 < nk; i2++) {
		if (sum) {
		    spec[i2][i1] += sf_cabsf(ctrace2[i2]);
		} else { 
		    spec[i2][i1] = sf_cabsf(ctrace2[i2])*scale;
		}
	    }
	}
#endif
	
	if (!sum) sf_floatwrite(spec[0],nw*nk,out);
    } /* i3 */

    if (sum) { /* normalize and output */
	scale /= n3;
	for (i2=0; i2 < nk; i2++) { 
	    for (i1=0; i1 < nw; i1++) {
		spec[i2][i1] *= scale;
	    }
	}
	sf_floatwrite(spec[0],nw*nk,out);
    }

    exit (0);
}
