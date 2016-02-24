/* Streaming PEF */
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

int main(int argc, char* argv[])
{
    bool inv, adj, linear, kn;
    int n1, n2, na, i1, i2, ia, *mask;
    float dd, da, dn, rn, eps;
    float *d, *a, *r, *d2, *r2;
    sf_file inp, pef, out, pat, known;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* inversion flag */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag (for linear operator) */

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);

    if (!sf_getint("na",&na)) sf_error("Need na=");
    /* PEF filter size (not including leading one) */
    if (na > n1) sf_error("Cannot handle na > n1");

    if (!sf_getfloat("eps",&eps)) sf_error("Need eps=");
    /* regularization */
    eps *= eps;

    d = sf_floatalloc(n1);
    r = sf_floatalloc(n1);
    a = sf_floatalloc(na);

    linear = (NULL != sf_getstring("pattern"));
    /* pattern data (for linear operator) */

    if (linear) {
	pat = sf_input("pattern");
	d2  = sf_floatalloc(n1);
	r2  = sf_floatalloc(n1);
    } else {
	pat = inp;
	d2 = NULL;
	r2 = NULL;
    }

    if (NULL != sf_getstring("pef")) {
	/* output PEF (optional) */
	pef = sf_output("pef");
	sf_putint(pef,"n1",na);
	sf_putint(pef,"n2",n1);
	sf_putint(pef,"n3",n2);
    } else {
	pef = NULL;
    }

    if (NULL != sf_getstring("known")) {
	/* known data locations (optional) */
	known = sf_input("known");
	if (SF_INT != sf_gettype(known)) sf_error("Need int type in known");
	mask = sf_intalloc(n1);
    } else {
	known = NULL;
    }

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(inv? r: d,n1,pat);
	if (NULL != known) sf_intread(mask,n1,known);

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

	dd = 0.0f;
	da = 0.0f;
	for (ia=0; ia < na; ia++) {
	    a[ia]=0.0f;
	    if (inv) {
		d[ia] = r[ia];
	    } else {
		r[ia] = d[ia];
	    }
	    dd += d[ia]*d[ia];
	}

	if (linear) {
	    for (ia=0; ia < na; ia++) {
		if (adj) {
		    d2[ia] = r2[ia];
		} else {
		    r2[ia] = d2[ia];
		}
	    }
	}

	if (NULL != pef) {
	    for (ia=0; ia < na; ia++) {
		sf_floatwrite(a,na,pef);
	    }
	}
	
	for (i1=na; i1 < n1; i1++) {
	    kn = (NULL == known) || (1==mask[i1]);

	    if (inv) {
		rn = kn? r[i1]/eps: 0.0f;
		dn = rn*(eps+dd)-da;
		d[i1] = dn;
	    } else {
		dn = kn? d[i1]: -da;
		rn = (dn+da)/(eps+dd);
		r[i1] = eps*rn;
	    }

	    for (ia=0; ia < na; ia++) {
		a[ia] -= rn*d[i1-1-ia];
	    }

	    if (NULL != pef) sf_floatwrite(a,na,pef);

	    dd += dn*dn - d[i1-na]*d[i1-na];	    
	    da = dn*a[0];
	    for (ia=1; ia < na; ia++) {
		da += a[ia]*d[i1-ia];
	    }

	    if (linear) {
		if (adj) {
		    d2[i1] += r2[i1];
		    for (ia=0; ia < na; ia++) {
			d2[i1-1-ia] += a[ia]*r2[i1];
		    }
		} else {
		    r2[i1] = d2[i1];
		    for (ia=0; ia < na; ia++) {
			r2[i1] += a[ia]*d2[i1-1-ia];
		    }
		}
	    }
	}

	if (linear) {
	    sf_floatwrite(adj? d2: r2,n1,out);
	} else {	
	    sf_floatwrite(inv? d: r,n1,out);
	}
    }

    exit(0);
}
