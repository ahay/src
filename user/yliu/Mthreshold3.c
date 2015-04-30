/* Automatic soft or hard thresholding. */
/*
  Copyright (C) 2015 Jiln University

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
    int i, nd, nr, ir;
    int n[SF_MAX_DIM], ndim;
    float *dat=NULL, *adat=NULL;
    float t, d, coef, median, median1, median2, sigma=0.;
    char *type, *dist;
    bool soft=false;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (NULL == (type=sf_getstring("type"))) type="soft";
    /* [soft,hard] thresholding type, the default is soft  */

    if (NULL == (dist=sf_getstring("dist"))) dist="gaussian";
    /* [gaussian,rayleigh] distribution type, the default is gaussian */

    switch(type[0]) {
	case 's':
	    sf_warning("Thresholding type=%c",type[0]);
	    soft = true;
	    break;
	case 'h':
	    sf_warning("Thresholding type=%c",type[0]);
	    soft = false;
	    break;
	default:
	    sf_error("Unknown thresholding type=%c",type[0]);
	    break;
    }

    switch(dist[0]) {
	case 'g':
	    sigma = 0.6745;
	    break;
	case 'r':
	    sigma = 0.4485;
	    break;
	default:
	    sf_error("Unknown distribution type=%c",dist[0]);
	    break;
    }

    ndim = sf_filedims(in,n);
    nd = 1;
    for (i=0; i < ndim-1; i++) {
	nd *= n[i];
    }

    nr = sf_leftsize(in,ndim-1);

    adat = sf_floatalloc(nd);
    dat = sf_floatalloc(nd);
    coef = sqrt(2*log10(nd));

    for (ir=0; ir < nr; ir++) {    
	sf_warning("slice %d of %d;", ir+1, nr);

	sf_floatread(dat,nd,in);
	for (i=0; i < nd; i++) {
	    adat[i] = dat[i];
	}

	/* compute automatic thresholding value */
	if (nd%2) { /* odd number */
	    median = sf_quantile((nd-1)/2,nd,adat);
	} else { /* even number */
	    median1 = sf_quantile(nd/2-1,nd,adat);
	    median2 = sf_quantile(nd/2,nd,adat);
	    median = 0.5*(median1+median2);
	}

	for (i=0; i < nd; i++) {
	    adat[i] = fabsf(dat[i]-median);
	}

	if (nd%2) { /* odd number */
	    median = sf_quantile((nd-1)/2,nd,adat);
	} else { /* even number */
	    median1 = sf_quantile(nd/2-1,nd,adat);
	    median2 = sf_quantile(nd/2,nd,adat);
	    median = 0.5*(median1+median2);
	}

	t = median * coef / sigma;

	for (i=0; i < nd; i++) {
	    d = dat[i];
	    if (d < -t) {
		if(soft) dat[i] = d+t;
	    } else if (d > t) {
		if(soft) dat[i] = d-t;
	    } else {
		dat[i] = 0.;
	    }
	}

	sf_floatwrite(dat,nd,out);
	
    }
    sf_warning(".");
    exit(0);
}

/* 	$Id$	 */
