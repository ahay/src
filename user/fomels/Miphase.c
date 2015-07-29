/* Smooth estimate of instantaneous frequency. */
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
#include <math.h>

#include <rsf.h>

#include "cdivn.h"

int main (int argc, char* argv[])
{
    int nh, n1,n2, i1,i2, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *trace, *hilb, *dtrace, *dhilb, *num, *den, *phase, a,b,c, mean, d1;
    sf_complex *cnum, *cden, *crat;
    char key[6];
    bool hertz, band, cmplx, verb;
    sf_file in, out;
	
    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
	
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    n1 = n[0];
    n12 = 1;
    for (i=0; i < dim; i++) {
		snprintf(key,6,"rect%d",i+1);
		if (!sf_getint(key,rect+i)) rect[i]=1;
		/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
		n12 *= n[i];
    }
    n2 = n12/n1;

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
	
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	
    trace = sf_floatalloc(n1);
    hilb = sf_floatalloc(n1);
    dtrace = sf_floatalloc(n1);
    dhilb = sf_floatalloc(n1);
	
    if (!sf_getbool("complex",&cmplx)) cmplx=false;
    /* if y, use complex-valued computations */
	
    if (cmplx) {
		cnum = sf_complexalloc(n12);
		cden = sf_complexalloc(n12);
		crat = sf_complexalloc(n12);
		num = den = NULL;
    } else {
		num = sf_floatalloc(n12);
		den = sf_floatalloc(n12);
		cnum = cden = crat = NULL;
    }
	
    phase = sf_floatalloc(n12);
	
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
	
    if (!sf_getint("order",&nh)) nh=100;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */
	
    if (!sf_getbool("hertz",&hertz)) hertz=false;
    /* if y, convert output to Hertz */
	
    if (!sf_getbool("band",&band)) band=false;
    /* if y, compute instantaneous bandwidth */
	
    sf_hilbert_init(n1, nh, c);
    sf_deriv_init(n1, nh, c);
	
    mean=0.;
    for (i=i2=0; i2 < n2; i2++) {
	if (verb) sf_warning("slice %d of %d;",i2+1,n2);

	sf_floatread(trace,n1,in);
	sf_hilbert(trace,hilb);
	
	if (band) {
	    for (i1=0; i1 < n1; i1++) {
		/* find envelope */
		trace[i1] = hypotf(trace[i1],hilb[i1]);
	    }
	    sf_deriv(trace,hilb);
	} else {
	    sf_deriv(trace,dtrace);
	    sf_deriv(hilb,dhilb);
	}
	
	if (cmplx) {
	    for (i1=0; i1 < nh; i1++, i++) {
		cnum[i] = sf_cmplx(0.,0.);
		cden[i] = sf_cmplx(0.,0.);
	    }	
	    
	    for (i1=nh; i1 < n1-nh; i1++, i++) {
		cnum[i] = sf_cmplx(dtrace[i1],dhilb[i1]);
		cden[i] = sf_cmplx( trace[i1], hilb[i1]);
		
		a = cabsf(cden[i]);
		mean += a*a;
	    }
	    
	    for (i1=n1-nh; i1 < n1; i1++, i++) {
		cnum[i] = sf_cmplx(0.,0.);
		cden[i] = sf_cmplx(0.,0.);
	    }
	} else {
	    for (i1=0; i1 < nh; i1++, i++) {
		num[i] = 0.;
		den[i] = 0.;
	    }	
	    
	    for (i1=nh; i1 < n1-nh; i1++, i++) {
		a = trace[i1];
		b = hilb[i1];
		if (band) {
		    num[i] = b;
		    den[i] = a;
		} else {
		    num[i] = a*dhilb[i1]-b*dtrace[i1];
		    den[i] = a*a+b*b;
		}
		mean += den[i]*den[i];
	    }
	    
	    for (i1=n1-nh; i1 < n1; i1++, i++) {
		num[i] = 0.;
		den[i] = 0.;
	    }
	} /* cmplx */
	
    } /* i2 */
    
    if (verb) sf_warning(".");
    
    mean = sqrtf(n12/mean);
    
    for (i=0; i < n12; i++) {
		if (cmplx) {
#ifdef SF_HAS_COMPLEX_H
			cnum[i] *= mean;
			cden[i] *= mean;  
#else
			cnum[i] = sf_crmul(cnum[i],mean);
			cden[i] = sf_crmul(cden[i],mean);
#endif
		} else {
			num[i] *= mean;
			den[i] *= mean;
		}
    }
	
    if (cmplx) {
		cdivn_init(dim, n12, n, rect, niter, true);
		cdivn (cnum, cden, crat);
		for (i=0; i < n12; i++) {
			phase[i] = cimagf(crat[i]);
		}
    } else {
		sf_divn_init(dim, n12, n, rect, niter, true);
		sf_divn (num, den, phase);
    }
	
    if (hertz) {
		/* convert to Hertz */    
		d1 = 1./(2.*SF_PI*d1);
		for (i=0; i < n12; i++) {
			phase[i] *= d1;
		}
    }
	
    sf_floatwrite(phase,n12,out);
	
    exit(0);
}

/* 	$Id: Menvelope.c 696 2004-07-06 23:17:31Z fomels $	 */
