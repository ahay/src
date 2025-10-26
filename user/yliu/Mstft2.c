/* Another version of Short-time Fourier transform (STFT) with overlap windows. */
/*
  Copyright (C) 2025 Jilin University

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
#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    bool inv, sym, opt, wind, verb;
    int n1, nt, nw, i, i1, i2, n2, j, ntw, overlap, m, step, nfft;
    float dw, ow, *p, *inp, d1, o1, wt, shift, k, sum;
    kiss_fft_cpx *pp, ce;
    sf_complex *outp;
    sf_file in = NULL, out = NULL;
    kiss_fftr_cfg cfg;
    
    sf_init(argc, argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("inv", &inv)) inv = false;
    /* if y, perform inverse transform */
    if (!sf_getbool("wind", &wind))	wind = false;
    /* if y, add Hanning window*/
    if (!sf_getint("ntw", &ntw)) ntw = 7;
    /* time-window length */
    if (!sf_getint("overlap", &overlap))
	overlap = ntw - 1;
    if (overlap >= ntw)
	sf_error("Need overlap length smaller than window length.");
    /* overlap length */
    
    if (!sf_getint("nfft", &nfft)) nfft = 0;
    /* FFT size */
    if (!sf_getbool("sym", &sym)) sym = false;
    /* if y, apply symmetric scaling to make the FFT operator Hermitian */
    if (!sf_getbool("opt", &opt)) opt = true;
    /* if y, determine optimal size for efficiency */
    
    if (!sf_getbool("verb", &verb))	verb = false;
    /* verbosity flag */
    
    if (ntw < 1) sf_error("Need positive integer input");
    m = ntw / 2;
    
    if (inv) {
	if (SF_COMPLEX != sf_gettype(in))
	    sf_error("Need complex input");
	sf_settype(out, SF_FLOAT);
    } else {
	if (SF_FLOAT != sf_gettype(in))
	    sf_error("Need float input");
	sf_settype(out, SF_COMPLEX);
    }
    
    if (!inv) {
	if (!sf_histint(in, "n1", &n1))
	    n1 = 1;
	if (!sf_histfloat(in, "d1", &d1))
	    d1 = 1.;
	if (!sf_histfloat(in, "o1", &o1))
	    o1 = 0.;
	n2 = sf_leftsize(in, 1);
	
	/* Time interval: (window step length) step
	   Output coefficients array (complex) Coef: [nt,nw]
	*/
	step = ntw - overlap;
	nt = (n1) / step;
	if (nt * step < n1)
	    nt += 1;
	/* determine wavenumber sampling (for real to complex FFT) */
	if (0 == nfft)
	    nfft = opt ? 2 * kiss_fft_next_fast_size((ntw + 1) / 2) : ntw;
	else
	    nfft = opt ? 2 * kiss_fft_next_fast_size((nfft + 1) / 2) : nfft;
	nw = nfft / 2 + 1;
	dw = 1. / (nfft * d1);
	
	if (verb) {
	    sf_warning("-----------------------------------------------.");
	    sf_warning("| Window Length | FFT Length | Overlap Length |.");
	    sf_warning("|%14d |%11d |%15d |.", ntw, nfft, overlap);
	    sf_warning("-----------------------------------------------.");
	}
	
	sf_shiftdim(in, out, 1);
	sf_putint(out, "n1", nt);
	sf_putfloat(out, "d1", step * d1);
	sf_putint(out, "n2", nw);
	sf_putfloat(out, "d2", dw);
	sf_putfloat(out, "o2", 0.);
	sf_putstring(out, "label2", "Frequency");
	sf_putstring(out, "unit2", "Hz");
    } else {
	if (!sf_histint(in, "n1", &n1))
	    n1 = 1;
	if (!sf_histfloat(in, "d1", &d1))
	    d1 = 1.;
	if (!sf_histfloat(in, "o1", &o1))
	    o1 = 0.;
	
	if (!sf_histint(in, "n2", &nw))
	    sf_error("No n2= in input");
	if (!sf_histfloat(in, "d2", &dw))
	    sf_error("No d2= in input");
	if (!sf_histfloat(in, "o2", &ow))
	    ow = 0.;
	
	step = ntw - overlap;
	n1 = n1 * step;
	nt = (n1) / step;
	if (0 == nfft)
	    nfft = opt ? 2 * kiss_fft_next_fast_size((ntw + 1) / 2) : ntw;
	else
	    nfft = opt ? 2 * kiss_fft_next_fast_size((nfft + 1) / 2) : nfft;
	d1 = d1 / step;
	n2 = sf_leftsize(in, 2);
	sf_unshiftdim(in, out, 2);
	sf_putint(out, "n1", n1);
	sf_putfloat(out, "d1", d1);
    }
    
    p = sf_floatalloc(nfft > ntw ? nfft : ntw);
    pp = (kiss_fft_cpx *)sf_complexalloc(nfft);
    cfg = kiss_fftr_alloc(nfft, inv ? 1 : 0, 0, 0);
    wt = sym ? 1. / sqrtf((float)ntw) : 1.0 / ntw;
    
    inp = sf_floatalloc(n1);
    outp = sf_complexalloc(nw * nt);
    
    for (i2 = 0; i2 < n2; i2++) { // loop over traces
	sf_warning("slice %d of %d;", i2 + 1, n2);
	if (!inv) {
	    sf_floatread(inp, n1, in);
	    for (i = 0; i < nt; i += 1)	{ // loop over subseries
		for (j = 0; j < ntw; j++) {
		    if (i * step + j - m < 0 || i * step + j - m >= n1) {
			p[j] = 0.;
		    } else {
			/* Hanning window: w(n) = 0.5*(1-cos(n/N)) */
			if (wind)
			    p[j] = inp[i * step + j - m] * 0.5 * (1. - cosf(2 * SF_PI * j / (ntw - 1))) * wt;
			else
			    p[j] = inp[i * step + j - m] * wt;
		    }
		}
		if (ntw < nfft)
		    for (i1 = ntw; i1 < nfft; i1++)
			p[i1] = 0.;
		
		kiss_fftr(cfg, p, pp);
		
		if (0. != o1) {
		    for (j = 0; j < nfft; j++) {
			shift = -2.0 * SF_PI * j * dw * o1;
			ce.r = cosf(shift);
			ce.i = sinf(shift);
			pp[j] = sf_cmul(pp[j], ce);
		    }
		}
		
		for (j = 0; j < nw; j++) {
		    outp[j * nt + i] = sf_cmplx(pp[j].r, pp[j].i);
		}
	    }
	    sf_complexwrite(outp, nt * nw, out);
	} else {
	    sf_complexread(outp, nt * nw, in);
	    
	    for (i = 0; i < n1; i++) {
		inp[i] = 0.;
	    }
	    
	    for (i = 0; i < nt; i++) {
		for (i1 = 0; i1 < nw; i1++) {
		    pp[i1].r = crealf(outp[i1 * nt + i]);
		    pp[i1].i = cimagf(outp[i1 * nt + i]);
		}
		for (i1 = nw; i1 < nfft; i1++) {
		    pp[i1].r = crealf(0.);
		    pp[i1].i = cimagf(0.);
		}
		if (0. != o1) {
		    for (i1 = 0; i1 < nw; i1++)	{
			shift = +2.0 * SF_PI * i1 * dw * o1;
			ce.r = cosf(shift);
			ce.i = sinf(shift);
			pp[i1] = sf_cmul(pp[i1], ce);
		    }
		}
		
		kiss_fftri(cfg, pp, p);
		
		for (i1 = 0; i1 < ntw; i1++) {
		    p[i1] *= wt;
		    if (wind) {
			p[i1] /= (0.5 - 0.5 * cosf(2 * SF_PI * i1 / (ntw - 1)));
		    }
		    if (i * step + i1 - m >= 0 && i * step + i1 - m < n1) {
			inp[i * step + i1 - m] += p[i1];
		    }
		}
	    }
	    
	    if (wind) {
		sum = 0.;
		for (i1 = 0; i1 < ntw; i1++) {
		    k = 0.;
		    k = (0.5 - 0.5 * cosf(2 * SF_PI * i1 / (ntw - 1)));
		    sum += k * k;
		}
		for (i = 0; i < n1; i++) {
		    inp[i] = inp[i] / sum;
		}
	    } else {
		for (i = 0; i < m; i++)	{
		    inp[i] = inp[i] / ((i + m + 1) * 1.);
		}
		for (i = m; i < n1 - m; i++) {
		    inp[i] = inp[i] / (ntw * 1.);
		}
		for (i = n1 - m; i < n1; i++) {
		    inp[i] = inp[i] / ((n1 - i + m) * 1.);
		}
	    }
	    
	    sf_floatwrite(inp, n1, out);
	}
    }
    sf_warning(".");
    
    free(pp);
    free(p);
    free(cfg);
    free(inp);
    free(outp);
    
    exit(0);
}

/* 	$Id$	 */
