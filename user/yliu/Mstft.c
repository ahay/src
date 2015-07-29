/* Short-time Fourier transform (STFT). */
/*
  Copyright (C) 2011 Jilin University
  
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

int main(int argc, char* argv[])
{
    bool inv, sym, opt;
    int n1, nt, nw, i, i1, i2, n2, j, ntw, m;
    float dw, ow, *p, *inp, d1, o1, wt, shift;
    kiss_fft_cpx *pp, ce;
    sf_complex *outp;
    sf_file in=NULL, out=NULL;
    kiss_fftr_cfg cfg;

    sf_init(argc, argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, perform inverse transform */
    if (!sf_getbool("sym",&sym)) sym=false;
    /* if y, apply symmetric scaling to make the FFT operator Hermitian */
    if (!sf_getbool("opt",&opt)) opt=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getint("ntw",&ntw)) ntw=7;
    /* time-window length */

    if (ntw < 1)  sf_error("Need positive integer input"); 
    if (ntw%2 == 0)  ntw = (ntw+1);
    m = (ntw-1)/2;

    if (inv) {
	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	sf_settype (out,SF_FLOAT);
    } else {
	if (SF_FLOAT   != sf_gettype(in)) sf_error("Need float input");
    	sf_settype (out,SF_COMPLEX);
    }

    if (!inv) {
	if (!sf_histint  (in,"n1",&n1)) n1=1;
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;
	n2 = sf_leftsize(in,1);
	
	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((ntw+1)/2): ntw;
	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);

	sf_shiftdim(in, out, 1);
	sf_putint(out,"n2",nw);
	sf_putfloat(out,"d2",dw);
	sf_putfloat(out,"o2",0.);
	sf_putstring(out,"label2","Frequency");
	sf_putstring(out,"unit2","Hz");

    } else {
	if (!sf_histint  (in,"n1",&n1)) n1=1;
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	if (!sf_histint  (in,"n2",&nw)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&ow)) ow=0.; 

	n2 = sf_leftsize(in,2);
	nt = 2*(nw-1);
	sf_unshiftdim(in, out, 2);
    }	

    p = sf_floatalloc(nt);
    pp = (kiss_fft_cpx*) sf_complexalloc(nw);
    cfg = kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
    wt = sym? 1./sqrtf((float) nt): 1.0/nt;

    inp = sf_floatalloc(n1);
    outp = sf_complexalloc(n1*nw);

    for (i2=0; i2 < n2; i2++)  {
	sf_warning("slice %d of %d;",i2+1,n2);
	if (!inv) {
	    sf_floatread (inp,n1,in);
	    for (i=0; i < n1; i++)  {
		for (j=0; j < ntw; j++) {
		    if (i+j-m < 0 || i+j-m >= n1) {
			p[j] = 0.;
		    } else {
			p[j] = inp[i+j-m];
		    }
		}
		
		if (sym) {
		    for (i1=0; i1 < ntw; i1++) {
			p[i1] *= wt;
		    }
		}
		
		for (i1=ntw; i1 < nt; i1++) {
		    p[i1]=0.0;
		}
		
		kiss_fftr (cfg,p,pp);
		
		if (0. != o1) {
		    for (i1=0; i1 < nw; i1++) {
			shift = -2.0*SF_PI*i1*dw*o1;
			ce.r = cosf(shift);
			ce.i = sinf(shift);
			pp[i1]=sf_cmul(pp[i1],ce);
		    }
		}
		for (i1=0; i1 < nw; i1++) {
		    outp[i1*n1+i] = sf_cmplx(pp[i1].r,pp[i1].i);
		}
	    }
	    sf_complexwrite(outp,n1*nw,out);
	} else {
	    sf_complexread(outp,n1*nw,in);

	    for (i=0; i < n1; i++)  {
		inp[i] = 0.;
	    }
	    for (i=0; i < n1; i++)  {
		for (i1=0; i1 < nw; i1++) {
		    pp[i1].r = crealf(outp[i1*n1+i]);
		    pp[i1].i = cimagf(outp[i1*n1+i]);
		}
		if (0. != o1) {
		    for (i1=0; i1 < nw; i1++) {
			shift = +2.0*SF_PI*i1*dw*ow;
			ce.r = cosf(shift);
			ce.i = sinf(shift);
			pp[i1]=sf_cmul(pp[i1],ce);
		    }
		}
		
		kiss_fftri(cfg,pp,p);
		
		for (i1=0; i1 < ntw; i1++) {
		    p[i1] *= wt;
		    if (i+i1-m >= 0 && i+i1-m < n1) {
			inp[i+i1-m] += p[i1];
		    }
		}
	    }

	    for (i=0; i < m; i++) {
		inp[i] = inp[i]/((i+m+1)*1.);
	    }
	    for (i=m; i < n1-m; i++)  {
		inp[i] = inp[i]/(ntw*1.);
	    }
	    for (i=n1-m; i < n1; i++) {
		inp[i] = inp[i]/((n1-i+m)*1.);
	    }
	    sf_floatwrite(inp,n1,out);

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

/* 	$Id: Mstft.c 7202 2011-05-03 02:13:48Z yang_liu $	 */
