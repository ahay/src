/* Missing data interpolation in 2-D using Fourier Bregman shaping iteration. */
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
    int i, niter, nfft, nw, nk, n1, n2, n12;
    int i1, i2, i3, n3, iter; 
    float dw, /* dk, */ d1, o1, d2, o2, wt, wk, shift;
    float *dd, *dd2=NULL, *d, *dd3=NULL, *m, perc,fact;
    float ddif;
    float *err=NULL, *refd=NULL, *temp=NULL;
    bool verb, *known, error;
    sf_complex *mt;
    kiss_fft_cpx **mm, ce, **fft, *ctrace, *ctrace2;
    kiss_fftr_cfg tfft, itfft;
    kiss_fft_cfg  xfft, ixfft;
    sf_file in, out, mask=NULL, res=NULL, ref=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&d1)) d1=0.004;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    /* determine frequency sampling (for real to complex FFT) */
    nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
    nw = nfft/2+1;
    dw = 1./(nfft*d1);

    /* determine wavenumber sampling (for complex FFT) */
    nk = kiss_fft_next_fast_size(n2*2);
    /* dk = 1./(2.*nk*d2); */
    /* k0 = -0.5/d2; */

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool("error",&error)) error = false;
    /* error verbosity flag */

    dd = sf_floatalloc(n12);
    d = sf_floatalloc(nfft);
    ctrace = (kiss_fft_cpx*) sf_complexalloc (nw);
    ctrace2 = (kiss_fft_cpx*) sf_complexalloc (nk);
    mm = (kiss_fft_cpx**) sf_complexalloc2(nw,nk);
    fft = (kiss_fft_cpx**) sf_complexalloc2(nw,nk);
    mt = sf_complexalloc(nw*nk);
    known = sf_boolalloc(n12);
    m = sf_floatalloc(n12);

    tfft = kiss_fftr_alloc(nfft,0,NULL,NULL);
    xfft = kiss_fft_alloc(nk,0,NULL,NULL);
    itfft = kiss_fftr_alloc(nfft,1,NULL,NULL);
    ixfft = kiss_fft_alloc(nk,1,NULL,NULL);
    ddif = 0.;

    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }

    if (error && (NULL != sf_getstring ("res")) && 
	(NULL != sf_getstring ("ref"))) {
	res = sf_output("res");
	ref = sf_input("ref");
	refd = sf_floatalloc(n12);
	temp = sf_floatalloc(n12);
    } else {
	res = NULL;
	ref = NULL;
    }

    if (!sf_getfloat("perc",&perc)) perc=99.;
    /* percentage for soft-thresholding */
    if (!sf_getfloat("fact",&fact)) fact=0.5;
    /* factor for soft-thresholding */
    
    sf_sharpen_init(nw*nk,perc,fact);
    dd2 = sf_floatalloc(n12);
    dd3 = sf_floatalloc(n12);
    
    wt = 1.0/nfft; /* FFT time scaling */ 
    wk = 1.0/nk;   /* FFT distance scaling */ 

    if (error && (NULL != sf_getstring ("res")) && 
	(NULL != sf_getstring ("ref"))) {
	err = sf_floatalloc(niter);
	sf_putint(res,"n1",niter);
	sf_putfloat(res,"d1",1);
	sf_putfloat(res,"o1",1);	
	sf_putint(res,"n2",1);
    }

    for (i3=0; i3 < n3; i3++) {
	if (verb) {
	    sf_warning("slice %d of %d",i3+1,n3);
	} else {
	    sf_warning("slice %d of %d;",i3+1,n3);
	}
	sf_floatread(dd,n12,in);
	if (error && (NULL != sf_getstring ("res")) && 
	    (NULL != sf_getstring ("ref"))) {
	    sf_floatread(refd,n12,ref);
	}

	if (NULL != mask) {
	    sf_floatread(m,n12,mask);
	    for (i=0; i < n12; i++) {
		known[i] = (bool) (m[i] != 0.);
	    }
	} else {
	    for (i=0; i < n12; i++) {
		known[i] = (bool) (dd[i] != 0.);
	    }
	}
	for (i1=0; i1 < n12; i1++) {
	    dd2[i1] = 0.;
	    dd3[i1] = 0.;
	}
	for (iter=0; iter < niter; iter++) { /* Outer iteration */
	    if (verb)
		sf_warning("Bregman iteration %d of %d;",iter+1,niter);
	    
	    for (i1=0; i1 < n12; i1++) {
		if (known[i1]) {
		    dd3[i1]= dd[i1]+dd3[i1]-dd2[i1]; 
		}
	    }

	    for (i1=0; i1 < n12; i1++) {
		if (known[i1]) {
		    dd2[i1]= 0.; 
		}
	    }
	    for (i1=0; i1 < n12; i1++) {
		dd2[i1]+= dd3[i1]; 
	    }

	    /* Forward 2-D FFT */
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    d[i1] = dd2[i2*n1+i1];
		}
		for (i1=n1; i1 < nfft; i1++) {
		    d[i1] = 0.;
		}		    
		kiss_fftr (tfft,d,ctrace);
		
		for (i1=0; i1 < nw; i1++) {
		    fft[i2][i1] = i2%2? 
			sf_cneg(ctrace[i1]): ctrace[i1];
		}
		if (0. != o1) {
		    for (i1=0; i1 < nw; i1++) {
			shift = -2.0*SF_PI*i1*dw*o1;
			ce.r = cosf(shift);
			ce.i = sinf(shift);
			fft[i2][i1]=sf_cmul(fft[i2][i1],ce);
		    }
		}
	    }
	    
	    for (i2=n2; i2 < nk; i2++) {
		for (i1=0; i1 < nw; i1++) {
		    fft[i2][i1].r = 0.;
		    fft[i2][i1].i = 0.;
		}
	    }
	    
	    for (i1=0; i1 < nw; i1++) {
		/* Fourier transform x to k */
		kiss_fft_stride(xfft,fft[0]+i1,ctrace2,nw);
		for (i2=0; i2 < nk; i2++) {
		    mm[i2][i1] = ctrace2[i2];
		}
	    } /* Forward 2-D FFT end */
	    for (i2=0; i2 < nk; i2++) {
		for (i1=0; i1 < nw; i1++) {
		    mt[i2*nw+i1] = sf_cmplx(mm[i2][i1].r,
					    mm[i2][i1].i);
		}
	    }
	    /* Thresholding */
	    sf_csharpen(mt);
	    sf_cweight_apply(nw*nk, mt);
	    
	    for (i2=0; i2 < nk; i2++) {
		for (i1=0; i1 < nw; i1++) {
		    mm[i2][i1].r = crealf(mt[i2*nw+i1]);
		    mm[i2][i1].i = cimagf(mt[i2*nw+i1]);
		}
	    }
	    
	    /* Inverse 2-D FFT */
	    for (i1=0; i1 < nw; i1++) {
		/* Fourier transform k to x */
		kiss_fft_stride(ixfft,mm[0]+i1,ctrace2,nw);
		
		for (i2=0; i2 < nk; i2++) {
		    fft[i2][i1] = sf_crmul(ctrace2[i2],
					   i2%2? -wk: wk);
		}
	    }
	    for (i2=0; i2 < n2; i2++) {
		if (0. != o1) {
		    for (i1=0; i1 < nw; i1++) {
			shift = +2.0*SF_PI*i1*dw*o1;
			ce.r = cosf(shift);
			ce.i = sinf(shift);
			fft[i2][i1]=sf_cmul(fft[i2][i1],ce);
		    }
		}
		kiss_fftri(itfft,fft[i2],d);
		for (i1=0; i1 < n1; i1++) {
		    dd2[i2*n1+i1] = d[i1]*wt;
		}
	    } /* Inverse 2-D FFT end */
	    
	    if ((error && (NULL != sf_getstring ("res")) && 
		 (NULL != sf_getstring ("ref")))) {
		for (i1=0; i1 < n12; i1++) {
		    temp[i1] = dd2[i1];
		}
		for (i1=0; i1 < n12; i1++) {
		    temp[i1] -= refd[i1];
		}
		
		ddif = cblas_snrm2(n12, temp, 1);
		err[iter] = sqrtf(ddif/n12);
	    }
	} /* End Linear Bregman interation */
	
	for (i1=0; i1 < n12; i1++) {
	    dd[i1] = dd2[i1];
	}
	
	sf_floatwrite (dd,n12,out);
	
	if (error && (NULL != sf_getstring ("res")) && 
	    (NULL != sf_getstring ("ref"))) {
	    sf_floatwrite (err,niter,res);
	}
    }
    
    sf_warning(".");
    
    exit(0);
}

/* 	$Id$	 */
