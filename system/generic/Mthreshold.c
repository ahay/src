/* Soft thresholding. 

November 2014 program of the month:
http://ahay.org/blog/2014/11/12/program-of-the-month-sfthreshold/
*/
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

int main(int argc, char* argv[])
{
    int i, n, n1;
    float *dat=NULL, *adat=NULL, t, pclip, d;
    sf_complex *cdat=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    n = sf_filesize(in);
    adat = sf_floatalloc(n);

    if (!sf_getfloat("pclip",&pclip)) sf_error("Need pclip=");
    /* percentage to clip */
    n1 = 0.5+n*(1.-0.01*pclip);
    if (n1 < 0) n1=0;
    if (n1 >= n) n1=n-1;

    if (SF_FLOAT == sf_gettype(in)) {
	dat = sf_floatalloc(n);
	sf_floatread(dat,n,in);
	for (i=0; i < n; i++) {
	    adat[i] = fabsf(dat[i]);
	}
    } else if (SF_COMPLEX == sf_gettype(in)) {
	cdat = sf_complexalloc(n);
	sf_complexread(cdat,n,in);
	for (i=0; i < n; i++) {
	    adat[i] = cabsf(cdat[i]);
	}
    } else {
	sf_error("Need float or complex input");
    }

    t = sf_quantile(n1,n,adat);

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

    exit(0);
}
