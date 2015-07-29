/* 2-D Soft thresholding. */
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

int main(int argc, char* argv[])
{
    int i, n, nc, n1, n2, n3, nr, ir;
    float *dat=NULL, *adat=NULL, t, pclip, d, thrd;
    bool verb;
    sf_complex *cdat=NULL;
    sf_file in=NULL, out=NULL, thr=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histint(in,"n3",&n3)) n3=1;
    n = n1*n2;
    nr = sf_leftsize(in,2);

    adat = sf_floatalloc(n);

    if (NULL != sf_getstring ("thr")) {
	thr = sf_input("thr");
	if (nr != sf_filesize(thr)) 
	    sf_error("Need the same number as leftsize(input,2)");
    } else {
	thr = NULL;
    }

    if (!sf_getfloat("pclip",&pclip)) pclip=99.;

    if (SF_FLOAT == sf_gettype(in)) {
	dat = sf_floatalloc(n);
    } else if (SF_COMPLEX == sf_gettype(in)) {
	cdat = sf_complexalloc(n);
    } else {
	sf_error("Need float or complex input");
    }

    for (ir=0; ir < nr; ir++) {    
	if (verb) sf_warning("slice %d of %d;", ir+1, nr);

	if (NULL != thr) {
	    sf_floatread(&thrd,1,thr);
	} else {
	    thrd = pclip;
	}
	/* percentage to clip */
	nc = 0.5+n*(1.-0.01*thrd);
	if (nc < 0) nc=0;
	if (nc >= n) nc=n-1;
	
	if (SF_FLOAT == sf_gettype(in)) {
	    sf_floatread(dat,n,in);
	    for (i=0; i < n; i++) {
		adat[i] = fabsf(dat[i]);
	    }
	} else if (SF_COMPLEX == sf_gettype(in)) {
	    sf_complexread(cdat,n,in);
	    for (i=0; i < n; i++) {
		adat[i] = cabsf(cdat[i]);
	    }
	} else {
	    sf_error("Need float or complex input");
	}
	
	t = sf_quantile(nc,n,adat);
	
	if (NULL != dat) {
	    for (i=0; i < n; i++) {
		d = dat[i];
		if (d < -t) {
		    dat[i] = d+t;
		} else if (d > t) {
		    dat[i] = d-t;
		} else {
		    dat[i] = 0.;
		}
	    }
	    sf_floatwrite(dat,n,out);
	} else {
	    for (i=0; i < n; i++) {
		d = cabsf(cdat[i]);
		if (d < -t) {
#ifdef SF_HAS_COMPLEX_H
		    cdat[i] *= (d+t)/d;
#else
		    cdat[i] = sf_crmul(cdat[i],(d+t)/d);
#endif
		} else if (d > t) {		
#ifdef SF_HAS_COMPLEX_H
		    cdat[i] *= (d-t)/d;
#else
		    cdat[i] = sf_crmul(cdat[i],(d-t)/d);
#endif
		} else {
		    cdat[i] = sf_cmplx(0.,0.);
		}
	    }
	    sf_complexwrite(cdat,n,out);
	}
    }
    sf_warning(".");
    exit(0);
}

/* 	$Id: Mthreshold2.c 6426 2010-07-22 00:28:23Z yang_liu $	 */
