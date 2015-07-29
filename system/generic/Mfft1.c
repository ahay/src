/* Fast Fourier Transform along the first axis. */
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

int main (int argc, char *argv[])
{
    bool inv, sym, opt;
    int n1, nt, nw, i1, i2, n2;
    float dw, *p, d1, o1, wt, shift;
    kiss_fft_cpx *pp, ce;
    char *label;
    sf_file in=NULL, out=NULL;

#ifdef SF_HAS_FFTW
    fftwf_plan cfg;
#else
    kiss_fftr_cfg cfg;
#endif

    sf_init(argc, argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, perform inverse transform */
    if (!sf_getbool("sym",&sym)) sym=false;
    /* if y, apply symmetric scaling to make the FFT operator Hermitian */
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */

    if (inv) {
	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	sf_settype (out,SF_FLOAT);
    } else {
	if (SF_FLOAT   != sf_gettype(in)) sf_error("Need float input");
    	sf_settype (out,SF_COMPLEX);
    }

    n2 = sf_leftsize(in,1);

    if (!inv) {
	if (!sf_histint  (in,"n1",&n1)) n1=1;
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);

	sf_putint  (out,"n1",nw);
	sf_putfloat(out,"o1",0.);
	sf_putfloat(out,"d1",dw);

	sf_putfloat(out,"fft_o1",o1);
	sf_putfloat(out,"fft_n1",n1);

	/* fix label */
	if (NULL != (label = sf_histstring(in,"label1"))) {
	    sf_putstring(out,"fft_label1",label);
	    if (!sf_fft_label(1,label,out))
		sf_putstring(out,"label1","Wavenumber");
	}
    } else {
	if (!sf_histint  (in,"n1",&nw)) sf_error("No n1= in input");
	if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"fft_o1",&o1)) o1=0.; 

	nt = 2*(nw-1);
	d1 = 1./(nt*dw);

	if (!opt || !sf_histint  (in,"fft_n1",&n1)) n1 = nt;

	sf_putint  (out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);

	/* fix label */
	if (NULL != (label = sf_histstring(in,"fft_label1"))) {
	    sf_putstring(out,"label1",label);
	} else if (NULL != (label = sf_histstring(in,"label1"))) {
	    (void) sf_fft_label(1,label,out);
	}
    }	
    sf_fft_unit(1,sf_histstring(in,"unit1"),out);

    p = sf_floatalloc(nt);
    pp = (kiss_fft_cpx*) sf_complexalloc(nw);

#ifdef SF_HAS_FFTW
    if (inv) {
	cfg = fftwf_plan_dft_c2r_1d(nt, (fftwf_complex *) pp, p,
				    FFTW_ESTIMATE);
    } else {
	cfg = fftwf_plan_dft_r2c_1d(nt, p, (fftwf_complex *) pp,
				    FFTW_ESTIMATE);
    }
    if (NULL == cfg) sf_error("FFTW failure.");
#else
    cfg = kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
#endif

    wt = sym? 1./sqrtf((float) nt): 1.0/nt;

    for (i2=0; i2 < n2; i2++) {
	if (!inv) {
	    sf_floatread (p,n1,in);

	    if (sym) {
		for (i1=0; i1 < n1; i1++) {
		    p[i1] *= wt;
		}
	    }

	    for (i1=n1; i1 < nt; i1++) {
		p[i1]=0.0;
	    }

#ifdef SF_HAS_FFTW
	    fftwf_execute(cfg);
#else
	    kiss_fftr (cfg,p,pp);
#endif

	    if (0. != o1) {
		for (i1=0; i1 < nw; i1++) {
		    shift = -2.0*SF_PI*i1*dw*o1;
		    ce.r = cosf(shift);
		    ce.i = sinf(shift);
		    pp[i1]=sf_cmul(pp[i1],ce);
		}
	    }

	    sf_floatwrite((float*) pp,2*nw,out);
	} else {
	    sf_floatread((float*) pp,2*nw,in);

	    if (0. != o1) {
		for (i1=0; i1 < nw; i1++) {
		    shift = +2.0*SF_PI*i1*dw*o1;
		    ce.r = cosf(shift);
		    ce.i = sinf(shift);
		    pp[i1]=sf_cmul(pp[i1],ce);
		}
	    }

#ifdef SF_HAS_FFTW
	    fftwf_execute(cfg);
#else
	    kiss_fftri(cfg,pp,p);
#endif

	    for (i1=0; i1 < n1; i1++) {
		p[i1] *= wt;
	    }

	    sf_floatwrite (p,n1,out);
	}
    }


    exit (0);
}

/* 	$Id: Mfft1.c 8881 2012-07-26 21:20:50Z sfomel $	 */
