/* Soft or hard thresholding using exact-value or percentile thresholding.
   When mode=soft and ifperc=y, sfthreshold1 is equal to sfthreshold.
*/
/*
  Copyright (C) 2013 University of Texas at Austin

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
    int i, n, n1, ifperc;
    float *dat=NULL, *adat=NULL, t, thr, d;
    char *type;
    sf_complex *cdat=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    n = sf_filesize(in);
    adat = sf_floatalloc(n);

    if (NULL == (type=sf_getstring("type"))) type="soft";
    /* [soft,hard] thresholding type, the default is soft  */

    if(!sf_getint("ifperc",&ifperc)) ifperc=1;
    /* 0, exact-value thresholding; 1, percentile thresholding. */

    if (!sf_getfloat("thr",&thr)) sf_error("Need thr=");
    /* thresholding level */
    
    if(ifperc==1)
    {
    n1 = 0.5+n*(1.-0.01*thr);
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
    }
    else 
	t=thr;

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
		if(type[0]=='s') cdat[i] *= (d+t)/d;
#else
		if(type[0]=='s') cdat[i] = sf_crmul(cdat[i],(d+t)/d);
#endif
	    } else if (d > t) {		
#ifdef SF_HAS_COMPLEX_H
		if(type[0]=='s') cdat[i] *= (d-t)/d;
#else
		if(type[0]=='s') cdat[i] = sf_crmul(cdat[i],(d-t)/d);
#endif
	    } else {
		cdat[i] = sf_cmplx(0.,0.);
	    }
	}
	sf_complexwrite(cdat,n,out);
    }

    exit(0);
}
