/* Adaptive multidimensional nonstationary PEF. */
/*
  Copyright (C) 2016 University of Texas at Austin
  
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

#include "bound.h"
#include "misinput.h"
#include "createhelix.h" 

int main(int argc, char* argv[])
{
    bool verb;
    int n[SF_MAX_DIM], n0[SF_MAX_DIM], rect[SF_MAX_DIM];
    int a[SF_MAX_DIM], center[SF_MAX_DIM], gap[SF_MAX_DIM];
    int ndim, dim, n123, n123s, i, ia, ns, i1, niter, na, i4, n4, *kk;
    float *d, *f, *dd;
    double mean;
    char *filename, key[6];
    sf_filter aa;
    sf_file in, filt, lag, res, mask;
 
    sf_init(argc,argv);

    in = sf_input("in");
    filt = sf_output("out");

    if (NULL == (filename = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */

    lag = sf_output(filename);
    sf_putstring(filt,"lag",filename);

    if (NULL != (filename = sf_getstring("res"))) {
	/* output residual (optional) */
	res = sf_output(filename);
    } else {
	res = NULL;
    }

    ndim = sf_filedims(in,n);

    if (!sf_getint("dim",&dim)) dim=ndim; /* number of dimensions */

    n4 = sf_leftsize(in,dim);
    
    sf_putints (lag,"n",n,dim);

    if (!sf_getints("a",a,dim)) sf_error("Need a=");

    if (!sf_getints("center",center,dim)) {
	for (i=0; i < dim; i++) {
	    center[i] = (i+1 < dim && a[i+1] > 1)? a[i]/2: 0;
	}
    }

    if (!sf_getint("na",&na)) na=0;
    /* filter size */

    if (0 == na) {
	if (!sf_getints("gap",gap,dim)) {
	    for (i=0; i < dim; i++) {
		gap[i] = 0;
	    }
	}
	
	aa = createhelix(dim, n, center, gap, a); /* allocate PEF */
	
	for (i=0; i < dim; i++) {	    
	    n0[i] = n[i];
	}
    } else {
	aa =  sf_allocatehelix (na);
	if (!sf_getints ("lags", aa->lag, na)) sf_error("Need lags=");
	if (!sf_getints ("n", n0, dim)) {
	    for (i=0; i < dim; i++) {	    
		n0[i] = n[i];
	    }
	}
    }

    n123 = 1;
    for (i=0; i < dim; i++) {
	n123 *= n[i];
	
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
    }

    dd = sf_floatalloc(n123);
    kk = sf_intalloc(n123);

    if (NULL != sf_getstring("maskin")) {
	/* optional input mask file */
	mask = sf_input("maskin");

	switch (sf_gettype(mask)) {
	    case SF_INT:
		sf_intread (kk,n123,mask);
		break;
	    case SF_FLOAT:
		sf_floatread (dd,n123,mask);
		for (i=0; i < n123; i++) {
		    kk[i] = (dd[i] != 0.0f);
		}
		break;
	    default:
		sf_error ("Wrong data type in maskin");
		break;
	}

	sf_fileclose (mask);
    } else {
	for (i=0; i < n123; i++) {
	    kk[i] = 1;
	}
    }

    bound (dim, n0, n, a, aa);
    find_mask(n123, kk, aa);

    if (NULL != sf_getstring("maskout")) {
	/* optional output mask file */
	mask = sf_output("maskout");

	for (i=0; i < n123; i++) {
	    kk[i] = aa->mis[i]? 0: 1;
	}
	
	sf_settype(mask,SF_INT);
	sf_intwrite (kk,n123,mask);
    }

    na = aa->nh;

    sf_settype(lag,SF_INT);
    sf_putint(lag,"n1",na);
    for (i=1; i < dim; i++) {
	sprintf(key,"n%d",i+1);
	sf_putint(lag,key,1);
    }
    sf_intwrite(aa->lag,na,lag);
    sf_fileclose(lag);

    snprintf(key,3,"n%d",dim+1);
    sf_putint(filt,key,na);
    sf_shiftdim(in, filt, dim+1);    
    
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    n123s = n123*na;
    
    d = sf_floatalloc(n123s);
    f = sf_floatalloc(n123s);
    
    sf_multidivn_init(na, dim, n123, n, rect, d, NULL, verb); 

    for (i4=0; i4 < n4; i4++) {

	sf_floatread(dd,n123,in);
	
	/* apply shifts: dd -> d */

	mean = 0.0f;
	for (i=ia=0; ia < na; ia++) {
	    ns = aa->lag[ia];
	    for (i1=0; i1 < n123; i1++,i++) {
		if (i1 < ns || aa->mis[i1]) {
		    d[i] = 0.0f;
		} else {
		    d[i] = dd[i1-ns];
		    mean += d[i]*d[i];
		}
	    }
	}

	if (mean == 0.0f) {
	    sf_floatwrite(d,n123s,filt);
	    continue;
	}
	
	mean = sqrt (n123s/mean);

	/* -> apply mask */
	
	for(i=0; i < n123s; i++) {
	    d[i] *= mean;
	}
	for(i1=0; i1 < n123; i1++) {
	    if (aa->mis[i1]) {
		dd[i1] = 0.0f;
	    } else {
		dd[i1] *= mean;
	    }
	}

	sf_multidivn (dd,f,niter);
	sf_floatwrite(f,n123s,filt);

	if (NULL != res) {
	    for (i=ia=0; ia < na; ia++) {
		for (i1=0; i1 < n123; i1++,i++) {
		    dd[i1] -= d[i]*f[i];
		}
	    }
	    for (i1=0; i1 < n123; i1++) {
		dd[i1] /= mean;
	    }
	    
	    sf_floatwrite(dd,n123,res);
	}	
    }
    
    exit(0);
}
