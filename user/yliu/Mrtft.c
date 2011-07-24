/* Ricker time-frequency transform. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include <math.h>

#include "nmultidivn.h"

int main(int argc, char* argv[])
{
    bool inv, verb;
    int i1, n1, iw, nt, nw, i2, n2, is, shift, rect0, niter, n12;
    int m[SF_MAX_DIM], *rect;
    float t, d1, w, w0, dw, mean, alpha;
    float *trace, *kbsc, *mkbsc=NULL,*sscc=NULL, *mm=NULL, *outp;
    sf_file in, out, mask, basis;
   
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (NULL != sf_getstring("basis")) {
	basis = sf_output("basis");
    } else {
	basis = NULL;
    }

    if (!inv) {
	if (!sf_getint("nw",&nw)) { /* number of frequencies */
	    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
	    nw = nt/2+1;
	    dw = 1./(nt*d1);
	    w0 = 0.;
	} else {
	    if (!sf_getfloat("dw",&dw)) {
		/* frequency step */
		nt = 2*kiss_fft_next_fast_size((n1+1)/2);
		dw = 1./(nt*d1);
	    }
	    if (!sf_getfloat("w0",&w0)) w0=0.;
	    /* first frequency */
	}
	n2 = sf_leftsize(in,1);
	sf_shiftdim(in, out, 2);
	sf_putint(out,"n2",nw);
	sf_putfloat(out,"d2",dw);
	sf_putfloat(out,"o2",w0);
	sf_putstring(out,"label2","Frequency");
	sf_putstring(out,"unit2","Hz");
	sf_putint(out,"n3",n1);
	sf_putfloat(out,"d3",d1);
	sf_putfloat(out,"o3",0.);
	sf_putstring(out,"label3","Shift");
	sf_putstring(out,"unit3","s");

	if (!sf_getint("rect",&rect0)) rect0=10;
	/* smoothing radius */
	if (!sf_getint("niter",&niter)) niter=100;
	/* number of inversion iterations */
	if (!sf_getfloat("alpha",&alpha)) alpha=0.;
	/* frequency adaptivity */

	for(i2=0; i2 < SF_MAX_DIM; i2 ++) {
	    m[i2] = 1;
	}
	m[0] = n1;
    } else {
	n2 = sf_leftsize(in,3);
	if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");
	if (!sf_histfloat(in,"o2",&w0)) sf_error("No o2= in input");
	sf_unshiftdim2(in, out, 2);
	sf_settype(out,SF_FLOAT);
    }

    if (NULL != basis) {
	sf_shiftdim(in, basis, 2);
	sf_putint(basis,"n2",nw);
	sf_putfloat(basis,"d2",dw);
	sf_putfloat(basis,"o2",w0);
	sf_putstring(basis,"label2","Frequency");
	sf_putstring(basis,"unit2","Hz");
	sf_putint(basis,"n3",n1);
	sf_putfloat(basis,"d3",d1);
	sf_putfloat(basis,"o3",0.);

    }
	
    n12 = n1*n1*nw;
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;

    trace = sf_floatalloc(n1);
    kbsc    = sf_floatalloc(n12);
    outp = sf_floatalloc(n1*nw);

    rect = sf_intalloc(n1*nw);
    for (i1=0; i1 < n1; i1++) {
	for (iw=0; iw < nw; iw++) {
	    rect[i1*nw+iw] = SF_MAX(1, (int) rect0/(1.0+alpha*iw/nw));
	}
    }

    if (!inv) {
	sscc = sf_floatalloc(n12);
	nmultidivn_init(n1*nw, 1, n1, m, rect, kbsc, (bool) (verb && (n2 < 500)));
    }
    
    if (NULL != sf_getstring("mask")) {
	mask = sf_input("mask");
	mm = sf_floatalloc(n1);
	mkbsc = sf_floatalloc(n12);
    } else {
	mask = NULL;
	mm = NULL;
	mkbsc = NULL;
    }

    /* sin and cos basis */
    for (is=0; is < n1; is++) {
	shift = is*d1;
	for (iw=0; iw < nw; iw++) {
	    w = w0 + iw*dw;
	    if (0.==w) { /* zero frequency */
		for (i1=0; i1 < n1; i1++) {
		    kbsc[is*n1*nw+iw*n1+i1] = 1.;
		}
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    t = i1*d1;
		    kbsc[is*n1*nw+iw*n1+i1] = (1.-2.*SF_PI*SF_PI*w*w*
					       (t-shift)*(t-shift))*
			expf(-SF_PI*SF_PI*w*w*(t-shift)*(t-shift));
		}
	    }
	}
    }
    
    if (NULL != mm) {
	for (iw=0; iw < n12; iw++) {
	    mkbsc[iw] = kbsc[iw];
	}
    }

    if (NULL != basis) {
	sf_floatwrite(kbsc,n12,basis);
    }

    mean = 0.;
    for (i1=0; i1 < n12; i1++) {
	mean += kbsc[i1]*kbsc[i1];
    }
    mean = sqrtf (mean/(n12));
    for (i1=0; i1 < n12; i1++) {
	kbsc[i1] /= mean;
    }
    
    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);
	if (NULL != mm) {
	    sf_floatread(mm,n1,mask);
	    for (is=0; is < n1; is++) {
		for (iw=0; iw < 2*nw; iw++) {
		    for (i1=0; i1 < n1; i1++) {
			kbsc[is*n1*nw+iw*n1+i1] = mkbsc[is*n1*nw+iw*n1+i1]
			    *mm[i1];
		    }
		}
	    }
	    mean = 0.;
	    for (i1=0; i1 < n12; i1++) {
		mean += kbsc[i1]*kbsc[i1];
	    }
	    mean = sqrtf (mean/(n12));
	    for (i1=0; i1 < n12; i1++) {
		kbsc[i1] /= mean;
	    }
	}

	if (!inv) {
	    sf_floatread(trace,n1,in);
	    if (NULL != mm) {
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] *= mm[i1];
		}
	    }
	    
	    for(i1=0; i1 < n1; i1++) {
		trace[i1] /= mean;
	    }
	    nmultidivn (trace,sscc,niter);
	    
	    for (i1=0; i1 < nw*n1; i1++) {
		outp[i1] = 0.;
	    }
	    for (is=0; is < n1; is++) {
		for (iw=0; iw < nw; iw++) {
		    for (i1=0; i1 < n1; i1++) {
			outp[iw*n1+i1] = sscc[is*n1*nw+iw*n1+i1];
		    }
		}
		sf_floatwrite(outp,n1*nw,out);
	    }
	} else {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] = 0.;
	    }
	    for (is=0; is < n1; is++) {
		sf_floatread(outp,n1*nw,in);
		for (iw=0; iw < nw; iw++) {
		    for (i1=0; i1 < n1; i1++) {
			trace[i1] += outp[iw*n1+i1]*kbsc[is*n1*nw+iw*n1+i1]
			    *mean;
			if (NULL != mm) trace[i1] *= mm[i1];
		    }
		}
	    }
	    sf_floatwrite(trace,n1,out);
	}
    }
    sf_warning(".");

    exit(0);
}

/* 	$Id$	 */
