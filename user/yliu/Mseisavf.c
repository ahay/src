/* 1-D amplitude versus frequency (AVF) analysis with 1-D seislet frames */
/*
  Copyright (C) 2009 University of Texas at Austin
   
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

#include "freqlets.h"

int main(int argc, char *argv[])
{
    int n1, i2, n2, iw, nw, nt, n1w, ncycle, niter;
    bool verb;
    float *w0, d1, scale, *thr;
    char *type;
    sf_complex *pp, *qq, *z0;
    sf_file in, out, w, t;
    int iter, i, i1, it, ip;
    sf_complex *q0, *p0, *p1, *q1;
    float qdif0=0., pdif0=0., qdif, pdif, pi;
    

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    w = sf_input("freq");
    t = sf_input("thr");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");

    if (!sf_histint(w,"n1",&nw)) sf_error("No n1= in freq");    
    if (SF_FLOAT == sf_gettype(w)) {
	w0 = sf_floatalloc(nw);
	z0 = NULL;
    } else if (SF_COMPLEX == sf_gettype(w)) {
	w0 = NULL;
	z0 = sf_complexalloc(nw);
    } else {
	sf_error("Need float or complex type in freq");
	w0 = NULL;
	z0 = NULL;
    }

    if (!sf_histint(t,"n1",&nt)) sf_error("No n1= in thre");    
    if (nt != nw) sf_error("number of threshold != number of frequency");
    thr = sf_floatalloc(nt);   

    n1w = n1*nw;
    scale = 1./nw;

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getint("ncycle",&ncycle)) ncycle=0;
    /* number of iterations */

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of Bregman iterations */
    
    n2 = sf_leftsize(in,1);
    sf_putint(out,"n2",nw);
    (void) sf_shiftdim(in, out, 2);

    pp = sf_complexalloc(n1);   /* data space */
    qq = sf_complexalloc(n1w);  /* model space */
    q0 = sf_complexalloc(n1w);
    q1 = sf_complexalloc(n1);
    p0 = sf_complexalloc(n1);
    p1 = sf_complexalloc(n1);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    /* sampling in the input file */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */
 
    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	qdif0=0.;
	pdif0=0.;
	if (NULL != w0) {
	    sf_floatread(w0,nw,w);
	} else {
	    sf_complexread(z0,nw,w);
	}
	sf_floatread(thr,nt,t);

	/* Forward Transform */
	freqlets_init(n1,d1,true,true,type[0],nw,w0,z0);
	sf_complexread(pp,n1,in);

	/* sf_csharpinv function + variable thresholding */
	for (i1=0; i1 < n1; i1++) {
	    p0[i1] = pp[i1];
	}
	
	for (iter=0; iter < niter; iter++) { /* outer iteration */
	    freqlets_lop(true,false,n1w,n1,qq,p0);
	    
	    for (i1=0; i1 < n1w; i1++) {
		q0[i1] = qq[i1];
	    }
	    for (i1=0; i1 < n1; i1++) {
		p1[i1] = p0[i1];
	    }
	    
	    for (i=0; i < ncycle; i++) { /* inner iteration */
		freqlets_lop(false,false,n1w,n1,qq,p1);	    
		for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    p1[i1] *= (-scale);
#else
		    p1[i1] = sf_crmul(pp[i1],-scale);
#endif
		}
		freqlets_lop(true,true,n1w,n1,qq,p1);
		
		for (i1=0; i1 < n1w; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    qq[i1] += q0[i1];
#else
		    qq[i1] = sf_cadd(qq[i1],q0[i1]);
#endif
		}

		/* Variable thresholding */
		for (it=0; it < nt; it++) {
		    sf_sharpen_init(n1,thr[it],0.5);
		    for (ip=0; ip < n1; ip++) {
			q1[ip] = qq[it*n1+ip];
		    }			
		    sf_csharpen(q1);
		    sf_cweight_apply(n1,q1);
		    for (ip=0; ip < n1; ip++) {
			qq[it*n1+ip] = q1[ip];
		    }			
		}
		
		if (verb) {		  	    
		    qdif = 0.;
		    for (i1=0; i1 < n1w; i1++) {
			qdif += cabsf(qq[i1]);
		    }
		    
		    if (0==i) {
			qdif0 = qdif;
			qdif=1.;
		    } else {
			qdif /= qdif0;
		    }
		    
		sf_warning("inner iteration %d mnorm: %f",i,qdif);
		}
	    } /* inner iteration */
	    
	    freqlets_lop(false,false,n1w,n1,qq,p1);
	    for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		p1[i1] *= (-scale);
#else
		p1[i1] = sf_crmul(pp[i1],-scale);
#endif	
	    }
	    
	    for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		p0[i1] += pp[i1] + p1[i1];
#else
		p0[i1] = sf_cadd(p0[i1],sf_cadd(pp[i1],p1[i1]));
#endif
	    }
	    
	    if (verb) {
		pdif = 0.;
		for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		    pi = cabsf(pp[i1]+p1[i1]);
#else
		    pi = cabsf(sf_cadd(pp[i1],p1[i1]));
#endif
		    pdif += pi*pi;
		}
		if (0==iter) {
		    pdif0 = pdif;
		    pdif=1.;
		} else {
		    pdif /= pdif0;
		}
		sf_warning("outer iteration %d dres: %f",iter,pdif);
	    }
	} /* outer iteration */
	
	/* Decompositon and Inverse Transform */
	freqlets_init(n1,d1,true,true,type[0],1,w0,z0);
	for (iw=0; iw < nw; iw++) {
	    if (NULL != w0) {
		freqlets_set(w0+iw,NULL);
	    } else {
		freqlets_set(NULL,z0+iw);
	    }
	    freqlets_lop(false,false,n1,n1,qq+iw*n1,pp);
	    sf_complexwrite(pp,n1,out);
	}
    }

    exit(0);
}

/* 	$Id$	 */
