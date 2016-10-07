/* Streaming PEF on a helix */
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
#include "createhelix.h"
#include "printfilter.h"

int main(int argc, char* argv[])
{
    bool inv, adj, linear;
    int n1, na, i1, ia, i, maxlag, jump, dim, lag0, ii[SF_MAX_DIM];
    int n[SF_MAX_DIM], n0[SF_MAX_DIM], center[SF_MAX_DIM], a[SF_MAX_DIM], gap[SF_MAX_DIM];
    float dd, da, dn, rn, eps;
    float *d, *r, *d2=NULL, *r2=NULL;
    sf_filter aa;
    sf_file inp, out, pat, lag;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* inversion flag */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag (for linear operator) */

    if (!sf_getint("jump",&jump)) jump=1;
    /* jump > 1 is used for trace interpolation */

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
	    na = aa->nh;

	    if (jump > 1) {
		lag0 = sf_cart2line(dim, a, center);
		for (ia=0; ia < na; ia++) {	/* sweep through the filter */
		    sf_line2cart(dim, a, ia+lag0+1, ii);
		    for (i=0; i < dim; i++) {
			ii[i] -= center[i];
		    }
		    ii[0] *= jump;  /* interlace on 1-axis */
		    aa->lag[ia] = sf_cart2line(dim, n0, ii);
		}
	    }
	} else {	  
	    aa =  sf_allocatehelix (na);
	    if (!sf_getints ("lags", aa->lag, na)) sf_error("Need lags=");	    	    
	}
    }

    a[0] *= jump;
    bound (dim, false, n0, n, a, aa); 
    a[0] /= jump;
    
    maxlag = 0;
    for (ia=0; ia < na; ia++) {
	if (aa->lag[ia] > maxlag) maxlag=aa->lag[ia];
    }

    if (!sf_getfloat("eps",&eps)) sf_error("Need eps=");
    /* regularization */
    eps *= eps;

    d = sf_floatalloc(n1);
    r = sf_floatalloc(n1);

    linear = (NULL != sf_getstring("pattern"));
    /* pattern data (for linear operator) */

    if (linear) {
	pat = sf_input("pattern");
	d2  = sf_floatalloc(n1);
	r2  = sf_floatalloc(n1);
    } else {
	pat = inp;
    }

    sf_floatread(inv? r: d,n1,pat);

    if (linear)  {
	if (adj) {
	    sf_floatread(r2,n1,inp);
	    for (i1=0; i1 < n1; i1++) {
		d2[i1] = 0.0f;
	    }
	} else {
	    sf_floatread(d2,n1,inp);
	}
    }
    
    for (ia=0; ia < na; ia++) {
	aa->flt[ia]=0.0f;
    }

    for (ia=0; ia < maxlag; ia++) {
	if (inv) {
	    d[ia] = r[ia];
	} else {
	    r[ia] = d[ia];
	}
    }
    
    if (linear) {
	for (ia=0; ia < maxlag; ia++) {
	    if (adj) {
		d2[ia] = r2[ia];
	    } else {
		r2[ia] = d2[ia];
	    }
	}
    }

    dd = 0.0f;
    da = 0.0f;    
    for (ia=0; ia < na; ia++) {
	dd += d[maxlag-aa->lag[ia]]*d[maxlag-aa->lag[ia]];
    }
    
    for (i1=maxlag; i1 < n1; i1++) {
	if (inv) {
	    rn = r[i1]/eps;
	    dn = rn*(eps+dd)-da;
	    d[i1] = dn;
	} else {
	    dn = d[i1];
	    rn = (dn+da)/(eps+dd);
	    r[i1] = eps*rn;
	}

	for (ia=0; ia < na; ia++) {
	    aa->flt[ia] -= rn*d[i1-aa->lag[ia]];
	}

	dd += dn*dn - d[i1-maxlag]*d[i1-maxlag];
	
	da = dn*aa->flt[0];
	for (ia=1; ia < na; ia++) {
	    da += aa->flt[ia]*d[i1+1-aa->lag[ia]];
	}

	if (linear) {
	    if (adj) {
		d2[i1] += r2[i1];
		for (ia=0; ia < na; ia++) {
		    d2[i1-aa->lag[ia]] += aa->flt[ia]*r2[i1];
		}
	    } else {
		r2[i1] = d2[i1];
		for (ia=0; ia < na; ia++) {
		    r2[i1] += aa->flt[ia]*d2[i1-aa->lag[ia]];
		}
	    }
	}
    }

    if (linear) {
	sf_floatwrite(adj? d2: r2,n1,out);
    } else {	
	sf_floatwrite(inv? d: r,n1,out);
    }

    exit(0);
}
