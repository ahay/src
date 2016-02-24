/* Missing data interpolation using streaming PEF */
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
    int n1, n2, na, i1, i2, ia, *mask;
    float dd, da, dn, rn, eps;
    float *d, *a; 
    sf_file inp, pef, out, known;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

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
    a = sf_floatalloc(na);

    if (NULL != sf_getstring("pef")) {
	/* output PEF (optional) */
	pef = sf_output("pef");
	sf_putint(pef,"n1",na);
	sf_putint(pef,"n2",n1);
	sf_putint(pef,"n3",n2);
    } else {
	pef = NULL;
    }
    /* known data locations (optional) */
    known = sf_input("known");

    if (SF_INT != sf_gettype(known)) sf_error("Need int type in known");
    mask = sf_intalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(d,n1,inp);
	sf_intread(mask,n1,known);

	dd = 0.0f;
	da = 0.0f;
	for (ia=0; ia < na; ia++) {
	    a[ia]=0.0f;
	    dd += d[ia]*d[ia];
	}

	if (NULL != pef) {
	    for (ia=0; ia < na; ia++) {
		sf_floatwrite(a,na,pef);
	    }
	}
	
	for (i1=na; i1 < n1; i1++) {
	    if (mask[i1]) {
		dn = d[i1];
		rn = (dn+da)/(eps+dd);
	    } else {
		dn = -da;
		rn = 0.0f;
		d[i1] = dn;
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
	}

	sf_floatwrite(d,n1,out);
    }
    
    exit(0);
}
