/* Streaming PEF with missing data */
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
    bool inv, adj, linear;
    int n1, n2, na, i1, i2, ia, neq, *m1, *m2;
    float dd, da, dn, rn, eps, *d2=NULL, *r2=NULL;
    float *d, *a, *r;
    sf_file inp, pef, out, minp, mout, pat;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    minp = sf_input("maskin");
    mout = sf_output("maskout");
    sf_settype(mout,SF_INT);
    
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

    m1 = sf_intalloc(n1);
    m2 = sf_intalloc(n1);
    
    if (NULL != sf_getstring("pef")) {
	/* output PEF (optional) */
	pef = sf_output("pef");
	sf_putint(pef,"n1",na);
	sf_putint(pef,"n2",n1);
	sf_putint(pef,"n3",n2);
    } else {
	pef = NULL;
    }

    linear = (NULL != sf_getstring("pattern"));
    /* pattern data (for linear operator) */

    if (linear) {
	pat = sf_input("pattern");
	d2  = sf_floatalloc(n1);
	r2  = sf_floatalloc(n1);
    } else {
	pat = inp;
    }

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(inv? r: d,n1,pat);
	sf_intread(m1,n1,minp);

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
	    m2[ia]=1;
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

	neq = na;
	
	for (i1=na; i1 < n1; i1++) {
	    if (inv) {
		rn = r[i1]/eps;
		dn = rn*(eps+dd)-da;
		d[i1] = dn;
	    } else {
		dn = d[i1];
		rn = (dn+da)/(eps+dd);
		r[i1] = eps*rn;
	    }

	    neq += m1[i1] - m1[i1-na];
	    
	    if (neq == na) {
		/* update PEF if enough equations */
		for (ia=0; ia < na; ia++) {
		    a[ia] -= rn*d[i1-1-ia];
		}
		m2[i1]=1;
	    } else {
		m2[i1]=0;
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

	sf_intwrite(m2,n1,mout);
    }

    exit(0);
}
