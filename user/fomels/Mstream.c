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
    int n1, n2, na, i1, i2, ia;
    float *d, *a, *r, dd, da, dn, rn, eps;
    sf_file data, pef, res;

    sf_init(argc,argv);
    data = sf_input("in");
    res = sf_output("out");
    pef = sf_output("pef");

    if (SF_FLOAT != sf_gettype(data)) sf_error("Need float input");
    if (!sf_histint(data,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(data,1);

    if (!sf_getint("na",&na)) sf_error("Need na=");
    /* PEF filter size (not including leading one) */
    if (na > n1) sf_error("Cannot handle na > n1");

    if (!sf_getfloat("eps",&eps)) sf_error("Need eps=");
    /* regularization */
    eps *= eps;

    d = sf_floatalloc(n1);
    r = sf_floatalloc(n1);
    a = sf_floatalloc(na);

    sf_putint(pef,"n1",na);
    sf_putint(pef,"n2",n1);
    sf_putint(pef,"n3",n2);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(d,n1,data);

	dd = 0.0f;
	for (ia=0; ia < na; ia++) {
	    a[ia]=0.0f;
	    dd += d[ia]*d[ia];
	}
	da = 0.0f;
	for (ia=0; ia < na; ia++) {
	    r[ia] = d[ia];
	    sf_floatwrite(a,na,pef);
	}
	
	for (i1=na; i1 < n1; i1++) {
	    rn = dn = d[i1];
	    for (ia=0; ia < na; ia++) {
		a[ia] -= (dn+da)/(eps+dd)*d[i1-1-ia];
		rn += a[ia]*d[i1-1-ia];
	    }
	    r[i1] = rn;
	    sf_floatwrite(a,na,pef);

	    dd += dn*dn - d[i1-na]*d[i1-na];
	    
	    da = dn*a[0];
	    for (ia=1; ia < na; ia++) {
		da += a[ia]*d[i1-ia];
	    }
	}

	sf_floatwrite(r,n1,res);
    }

    exit(0);
}
