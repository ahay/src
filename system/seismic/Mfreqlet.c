/* 1-D seislet frame */
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
#include "freqlets.h"

static void datawrite(int n1, float scale, sf_complex *pp, sf_file out);

int main(int argc, char *argv[])
{
    int n1, i2, n2, iw, nw, n1w, ncycle, niter;
    bool inv, verb, decomp;
    float *w0=NULL, d1, perc, scale;
    char *type=NULL;
    sf_complex *pp=NULL, *qq=NULL, *z0=NULL;
    sf_file in=NULL, out=NULL, w=NULL;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    w = sf_input("freq");

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
    n1w = n1*nw;
    scale = 1./nw;

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    if (!sf_getbool("decomp",&decomp)) decomp=false;
    /* do decomposition */

    if (!sf_getint("ncycle",&ncycle)) ncycle=0;
    /* number of iterations */

    if (!sf_getint("niter",&niter)) niter=1;
    /* number of Bregman iterations */

    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    if (inv) {
	n2 = sf_leftsize(in,2);
	if (!decomp) sf_unshiftdim(in, out, 2);	
    } else {
	n2 = sf_leftsize(in,1);
	sf_putint(out,"n2",nw);
	(void) sf_shiftdim(in, out, 2);
    }

    pp = sf_complexalloc(n1);   /* data space */
    qq = sf_complexalloc(n1w);  /* model space */

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    /* sampling in the input file */

    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    freqlets_init(n1,d1,true,true,type[0],decomp? 1: nw,w0,z0);

    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	if (NULL != w0) {
	    sf_floatread(w0,nw,w);
	} else {
	    sf_complexread(z0,nw,w);
	}

	if (inv) {
	    sf_complexread(qq,n1w,in);

	    if (decomp) {
		for (iw=0; iw < nw; iw++) {
		    if (NULL != w0) {
			freqlets_set(w0+iw,NULL);
		    } else {
			freqlets_set(NULL,z0+iw);
		    }
		    freqlets_lop(false,false,n1,n1,qq+iw*n1,pp);
		    datawrite(n1,1.,pp,out);
		}
	    } else {
		freqlets_lop(false,false,n1w,n1,qq,pp);
		datawrite(n1,scale,pp,out);
	    }
	} else {
	    sf_complexread(pp,n1,in);
	    sf_csharpinv(freqlets_lop,
			 scale,niter,ncycle,perc,verb,n1w,n1,qq,pp,true);
	    sf_complexwrite(qq,n1w,out);
	}
    }

    exit(0);
}


static void datawrite(int n1, float scale, sf_complex *pp, sf_file out) 
{
    int i1;

    for (i1=0; i1 < n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	pp[i1] *= scale;
#else
	pp[i1] = sf_crmul(pp[i1],scale);
#endif
    }
    sf_complexwrite(pp,n1,out);
}
