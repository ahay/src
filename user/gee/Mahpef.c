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
    int ndim, dim, n123, n123s, i, niter, na, i4, n4, *kk;
    float *d, *f, *dd, mean;
    char *lagfile, key[6];
    sf_filter aa;
    sf_file in, filt, lag;
 
    sf_init(argc,argv);

    in = sf_input("in");
    filt = sf_output("out");

    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */

    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(filt,"lag",lagfile);

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

    for (i=0; i < n123; i++) {
	kk[i] = 1;
    }

    bound (dim, n0, n, a, aa);
    find_mask(n123, kk, aa);

    na = aa->nh;

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
	/* apply mask */
	
	mean = 0.;
	for(i=0; i < n123s; i++) {
	    mean += d[i]*d[i];
	}
	if (mean == 0.) {
	    sf_floatwrite(d,n123s,filt);
	    continue;
	}
	
	mean = sqrtf (mean/n123s);
	
	for(i=0; i < n123s; i++) {
	    d[i] /= mean;
	}
	for(i=0; i < n123; i++) {
	    dd[i] /= mean;
	}

	sf_multidivn (dd,f,niter);
	sf_floatwrite(f,n123s,filt);
	
    }
    
    exit(0);
}
