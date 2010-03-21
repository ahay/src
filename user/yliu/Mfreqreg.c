/* Local frequency interpolation. */
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

int main(int argc, char* argv[])
{
    int i1, n1, iw, nt, nw, i2, n2, n12;
    int rect, niter, m[SF_MAX_DIM], rec[SF_MAX_DIM];
    float t, d1, w, w0, dw, *trace, *bsc, *kbsc, *sscc, *mm, *reg, mean;
    bool verb;
    sf_filter aa=NULL;
    sf_file in, out, mask;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    n2 = sf_leftsize(in,1);
    
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
    
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;
    n12 = 2*n1*nw;
    
    trace = sf_floatalloc(n1);
    reg = sf_floatalloc(n1);
    bsc    = sf_floatalloc(n12);
    kbsc    = sf_floatalloc(n12);
    sscc    = sf_floatalloc(n12);
    
    if (!sf_getint("rect",&rect)) rect=10;
    /* smoothing radius */
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of inversion iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    
    for(i2=0; i2 < SF_MAX_DIM; i2 ++) {
	m[i2] = 1;
	rec[i2] = 1;
    }
    m[0] = n1;
    rec[0] = rect;
    
    sf_multidivn_init(2*nw, 1, n1, m, rec, kbsc, aa, false); 
    
    if (NULL != sf_getstring("mask")) {
	mask = sf_input("mask");
	mm = sf_floatalloc(n1);
    } else {
	mask = NULL;
	mm = NULL;
    }
    
    if (NULL != mm) {
        sf_floatread(mm,n1,mask);
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] *= mm[i1];
	}
    }
    /* sin and cos basis */
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    t = i1*d1;
	    kbsc[iw*n1+i1] = sinf(w*t);
	    if (NULL != mm) {
	        kbsc[iw*n1+i1] *= mm[i1];
	    }
	}
    }
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    t = i1*d1;
	    kbsc[(iw+nw)*n1+i1] = cosf(w*t);
	    if (NULL != mm) {
	        kbsc[(iw+nw)*n1+i1] *= mm[i1];
	    }
	}
    }
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    t = i1*d1;
	    bsc[iw*n1+i1] = sinf(w*t);
	}
    }
    for (iw=0; iw < nw; iw++) {
        w = w0 + iw*dw;
	for (i1=0; i1 < n1; i1++) {
	    t = i1*d1;
	    bsc[(iw+nw)*n1+i1] = cosf(w*t);
	}
    }
    mean = 0.;
    for (i2=0; i2 < n12; i2++) {
        mean += kbsc[i2]*kbsc[i2];
    }
    mean = sqrtf (mean/(n12));
    for (i1=0; i1 < n12; i1++) {
        kbsc[i1] /= mean;
    }
    
    for (i2=0; i2 < n2; i2++) {
	if (verb) sf_warning("slice %d of %d",i2+1,n2);
	sf_floatread(trace,n1,in);
	for (i1=0; i1 < n1; i1++) {
	    reg[i1] = 0.;
	}
	
	for(i1=0; i1 < n1; i1++) {
	    trace[i1] /= mean;
	}
	sf_multidivn (trace,sscc,niter);
	
	for (iw=0; iw < nw; iw++) {
	    for (i1=0; i1 < n1; i1++) {
	        reg[i1] +=sscc[iw*n1+i1]*bsc[iw*n1+i1];
	    }
	}
	for (iw=0; iw < nw; iw++) {
	    for (i1=0; i1 < n1; i1++) {
	        reg[i1] +=sscc[(iw+nw)*n1+i1]*bsc[(iw+nw)*n1+i1];
	    }
	}
	sf_floatwrite(reg,n1,out);
    }
    exit(0);
}

/* 	$Id$	 */
