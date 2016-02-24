/* Missing data interpolating using streaming PEF on a helix */
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
#include <math.h>
#include <time.h>

#include <rsf.h>

#include "bound.h"	
#include "createhelix.h"
#include "printfilter.h"

int main(int argc, char* argv[])
{
    int n1, na, i1, ia, i, maxlag, seed;
    int dim, *mask;
    int n[SF_MAX_DIM], n0[SF_MAX_DIM], center[SF_MAX_DIM], a[SF_MAX_DIM], gap[SF_MAX_DIM];
    float dd, da, dn, rn, eps, var;
    float *d;
    sf_filter aa;
    sf_file inp, out, lag, known;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");

    dim = sf_filedims(inp,n);
    n1 = 1;
    for (i=0; i < dim; i++) {
	n1 *= n[i];
    }

    if (NULL != sf_getstring("lag")) {
	lag = sf_input("lag");
	if (SF_INT != sf_gettype(lag)) sf_error("Need int lag");
	sf_histint(lag,"n1",&na);
	sf_histints(lag,"n",n0,dim);

	aa =  sf_allocatehelix (na);
	sf_intread(aa->lag,na,lag);
    } else {
	if (!sf_getint("na",&na)) na=0;
	/* PEF filter size (not including leading one) */
	if (na > n1) sf_error("Cannot handle na > n1");
	
	if (!sf_getints("a",a,dim)) sf_error("Need a=");
	/* filter shape */
	
	if (!sf_getints ("n", n0, dim)) {
	    for (i=0; i < dim; i++) {	    
		n0[i] = n[i];
	    }
	}
	
	if (0 == na) {
	    for (i=0; i < dim; i++) {
		center[i] = (i+1 < dim && a[i+1] > 1)? a[i]/2: 0;
		gap[i] = 0;
	    }
	    aa = createhelix(dim, n0, center, gap, a); /* allocate PEF */
	    print(dim, n0, center, a, aa);             /* print filter */
	    na = aa->nh;
	} else {	  
	    aa =  sf_allocatehelix (na);
	    if (!sf_getints ("lags", aa->lag, na)) sf_error("Need lags=");	    	    
	}
    }
	
    bound (dim, n0, n, a, aa); 
	
    maxlag = 0;
    for (ia=0; ia < na; ia++) {
	if (aa->lag[ia] > maxlag) maxlag=aa->lag[ia];
    }

    if (!sf_getfloat("eps",&eps)) sf_error("Need eps=");
    /* regularization */
    eps *= eps;

    d = sf_floatalloc(n1);

    sf_floatread(d,n1,inp);

    /* known data locations (optional) */
    known = sf_input("known");

    if (SF_INT != sf_gettype(known)) sf_error("Need int type in known");

    mask = sf_intalloc(n1);
    sf_intread(mask,n1,known);
    
    for (ia=0; ia < na; ia++) {
	aa->flt[ia]=0.0f;
    }

    if (!sf_getfloat("var",&var)) var=0.0f;
    /* noise variance */
    var = sqrtf(var);

    if (!sf_getint("seed",&seed)) seed = time(NULL);
    /* random seed */
    init_genrand((unsigned long) seed);

    dd = 0.0f;
    da = 0.0f;    
    for (ia=0; ia < na; ia++) {
	dd += d[maxlag-aa->lag[ia]]*d[maxlag-aa->lag[ia]];
    }
    
    for (i1=maxlag; i1 < n1; i1++) {
	if (mask[i1]) {
	    dn = d[i1];
	    rn = (dn+da)/(eps+dd);
	} else {
	    rn = var*sf_randn_one_bm()/eps;
	    dn = rn*(eps+dd)-da;
	    d[i1] = dn;
	}

	for (ia=0; ia < na; ia++) {
	    aa->flt[ia] -= rn*d[i1-aa->lag[ia]];
	}

	dd += dn*dn - d[i1-maxlag]*d[i1-maxlag];
	
	da = dn*aa->flt[0];
	for (ia=1; ia < na; ia++) {
	    da += aa->flt[ia]*d[i1+1-aa->lag[ia]];
	}
    }

    sf_floatwrite(d,n1,out);

    exit(0);
}
